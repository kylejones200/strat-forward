//! OpenLineage. Events produced by real CLI runs are validated against the official OpenLineage
//! 2-0-2 schema and standard facet schemas vendored in tests/fixtures/openlineage, plus the custom facet
//! schema in schemas/openlineage. Every facet is validated against the schema its own `_schemaURL` names.
//! The HTTP transport is checked with an in-process responder.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpListener;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::{Arc, Mutex};

use jsonschema::{Draft, Resource};
use strat_forward::Config;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const RUN_ID: &str = "7a1d2b3c-4e5f-4a6b-8c7d-9e0f1a2b3c4d";

fn schemas() -> BTreeMap<String, Value> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut files: Vec<PathBuf> = std::fs::read_dir(root.join("tests/fixtures/openlineage"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "json"))
        .collect();
    files.push(root.join("schemas/openlineage/StratigraphyFacets.json"));
    files
        .into_iter()
        .map(|p| {
            let v: Value = serde_json::from_slice(&std::fs::read(&p).unwrap()).unwrap();
            (v["$id"].as_str().unwrap().to_string(), v)
        })
        .collect()
}

fn validate(reference: &str, instance: &Value) -> Vec<String> {
    let mut opts = jsonschema::options()
        .with_draft(Draft::Draft202012)
        .should_validate_formats(true);
    for (id, schema) in schemas() {
        opts = opts.with_resource(
            id,
            Resource::from_contents(schema).expect("schema resource"),
        );
    }
    let validator = opts
        .build(&json!({ "$ref": reference }))
        .unwrap_or_else(|e| panic!("cannot resolve {reference}: {e}"));
    validator
        .iter_errors(instance)
        .map(|e| format!("{reference} at {}: {e}", e.instance_path))
        .collect()
}

/// Validates the event, then every facet against its own `_schemaURL`.
fn assert_valid(event: &Value) {
    let mut errors = validate(
        "https://openlineage.io/spec/2-0-2/OpenLineage.json#/$defs/RunEvent",
        event,
    );
    let mut facets = Vec::new();
    let collect = |v: &Value, out: &mut Vec<Value>| {
        if let Some(map) = v.as_object() {
            out.extend(map.values().cloned());
        }
    };
    collect(&event["run"]["facets"], &mut facets);
    collect(&event["job"]["facets"], &mut facets);
    for ds in event["inputs"]
        .as_array()
        .unwrap()
        .iter()
        .chain(event["outputs"].as_array().unwrap())
    {
        collect(&ds["facets"], &mut facets);
        collect(&ds["outputFacets"], &mut facets);
    }
    assert!(!facets.is_empty());
    for f in &facets {
        let url = f["_schemaURL"].as_str().expect("facet without _schemaURL");
        errors.extend(validate(url, f));
    }
    assert!(
        errors.is_empty(),
        "invalid OpenLineage:\n{}",
        errors.join("\n")
    );
}

fn small_config(dir: &Path, nx: usize) -> PathBuf {
    let mut c = Config::example();
    c.length_m = 10_000.0;
    c.nx = nx;
    c.duration_yr = 40_000.0;
    c.record_dt_yr = 10_000.0;
    c.supply.flux_m2_per_yr = 1.0;
    c.initial.basement = strat_forward::config::Basement::Ramp {
        top_m: 5.0,
        bottom_m: -60.0,
    };
    let p = dir.join("run.json");
    std::fs::write(&p, serde_json::to_string_pretty(&c).unwrap()).unwrap();
    p
}

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("openlineage-{name}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn cli(config: &Path, out: &Path, env: &[(&str, &str)]) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_strat-forward"));
    cmd.arg("run")
        .arg("--config")
        .arg(config)
        .arg("--out")
        .arg(out)
        .args(["--run-id", RUN_ID]);
    cmd.env_remove("OPENLINEAGE_URL")
        .env_remove("OPENLINEAGE_API_KEY");
    for (k, v) in env {
        cmd.env(k, v);
    }
    cmd.output().unwrap()
}

fn events(out: &Path) -> Vec<Value> {
    serde_json::from_slice(&std::fs::read(out.join("openlineage-events.json")).unwrap()).unwrap()
}

fn sha256(path: &Path) -> String {
    Sha256::digest(std::fs::read(path).unwrap())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

#[test]
fn a_completed_run_emits_valid_start_and_complete_events() {
    let dir = workdir("complete");
    let config = small_config(&dir, 20);
    let out = dir.join("out");
    let o = cli(&config, &out, &[]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let ev = events(&out);
    assert_eq!(ev.len(), 2);
    for e in &ev {
        assert_valid(e);
        assert_eq!(e["run"]["runId"], RUN_ID);
        assert_eq!(e["job"]["name"], "strat-forward.run");
    }
    let (start, done) = (&ev[0], &ev[1]);
    assert_eq!(start["eventType"], "START");
    assert_eq!(done["eventType"], "COMPLETE");
    assert!(start["outputs"].as_array().unwrap().is_empty());

    // The input is the config file, versioned by content.
    let input = &done["inputs"][0];
    assert_eq!(input["namespace"], "file");
    assert_eq!(input["name"], config.display().to_string());
    assert_eq!(
        input["facets"]["version"]["datasetVersion"],
        format!("sha256:{}", sha256(&config))
    );

    // Every output names a file that exists, with its content hash and size.
    let outputs = done["outputs"].as_array().unwrap();
    let names: Vec<&str> = outputs
        .iter()
        .map(|d| d["name"].as_str().unwrap())
        .collect();
    for want in [
        "summary.json",
        "ledger.json",
        "stratigraphy.json",
        "prov.jsonld",
    ] {
        assert!(
            names.iter().any(|n| n.ends_with(want)),
            "missing output {want}"
        );
    }
    for d in outputs {
        let p = Path::new(d["name"].as_str().unwrap());
        assert_eq!(
            d["facets"]["version"]["datasetVersion"],
            format!("sha256:{}", sha256(p))
        );
        assert_eq!(
            d["outputFacets"]["outputStatistics"]["size"],
            std::fs::metadata(p).unwrap().len()
        );
    }

    let facet = &done["run"]["facets"]["stratigraphy_run"];
    assert_eq!(facet["config_source"], "file");
    assert!(facet["max_relative_mass_error"].as_f64().unwrap() < 1e-10);
    assert!(facet["ledger_replay_max_diff_m"].as_f64().unwrap() < 1e-9);

    // The PROV document describes the same run.
    let prov: Value =
        serde_json::from_slice(&std::fs::read(out.join("prov.jsonld")).unwrap()).unwrap();
    let run_iri = format!("urn:uuid:{RUN_ID}");
    assert!(prov["@graph"]
        .as_array()
        .unwrap()
        .iter()
        .any(|n| n["@id"] == run_iri.as_str()));
}

#[test]
fn a_failed_run_emits_a_valid_fail_event() {
    let dir = workdir("fail");
    let config = small_config(&dir, 2); // nx < 3 is rejected by validation
    let out = dir.join("out");
    let o = cli(&config, &out, &[]);
    assert_eq!(o.status.code(), Some(1));
    // The target is never created by a failed run; its files and events are kept beside it.
    assert!(
        !out.exists(),
        "a failed run must not create the output directory"
    );
    let ev = events(&dir.join(format!("out.failed-{RUN_ID}")));
    assert_eq!(ev.len(), 2);
    ev.iter().for_each(assert_valid);
    assert_eq!(ev[1]["eventType"], "FAIL");
    let msg = ev[1]["run"]["facets"]["errorMessage"]["message"]
        .as_str()
        .unwrap();
    assert!(msg.contains("nx"), "{msg}");
}

/// Requests a test server saw: request line, headers, body.
type Requests = Arc<Mutex<Vec<(String, Vec<(String, String)>, Vec<u8>)>>>;

/// Accepts `n` requests, answering each with `status`, and records request line, headers and body.
fn responder(n: usize, status: u16) -> (String, Requests) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    let seen = Arc::new(Mutex::new(Vec::new()));
    let log = seen.clone();
    std::thread::spawn(move || {
        for _ in 0..n {
            let Ok((stream, _)) = listener.accept() else {
                return;
            };
            let mut reader = BufReader::new(stream.try_clone().unwrap());
            let mut request_line = String::new();
            reader.read_line(&mut request_line).unwrap();
            let mut headers = Vec::new();
            loop {
                let mut line = String::new();
                reader.read_line(&mut line).unwrap();
                let line = line.trim_end();
                if line.is_empty() {
                    break;
                }
                let (k, v) = line.split_once(':').unwrap();
                headers.push((k.trim().to_ascii_lowercase(), v.trim().to_string()));
            }
            let len: usize = headers
                .iter()
                .find(|(k, _)| k == "content-length")
                .map(|(_, v)| v.parse().unwrap())
                .unwrap_or(0);
            let mut body = vec![0; len];
            reader.read_exact(&mut body).unwrap();
            log.lock()
                .unwrap()
                .push((request_line.trim_end().to_string(), headers, body));
            let mut out = stream;
            write!(
                out,
                "HTTP/1.1 {status} X\r\ncontent-length: 0\r\nconnection: close\r\n\r\n"
            )
            .unwrap();
        }
    });
    (url, seen)
}

#[test]
fn events_are_posted_with_the_openlineage_http_transport() {
    let dir = workdir("http");
    let config = small_config(&dir, 20);
    let out = dir.join("out");
    let (url, seen) = responder(2, 200);
    let o = cli(
        &config,
        &out,
        &[
            ("OPENLINEAGE_URL", &url),
            ("OPENLINEAGE_API_KEY", "test-key"),
        ],
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let file = events(&out);
    let got = seen.lock().unwrap();
    assert_eq!(got.len(), 2);
    for ((line, headers, body), expected) in got.iter().zip(&file) {
        assert_eq!(line, "POST /api/v1/lineage HTTP/1.1");
        assert!(headers.contains(&("authorization".into(), "Bearer test-key".into())));
        let posted: Value = serde_json::from_slice(body).unwrap();
        assert_eq!(&posted, expected, "the posted event is the recorded one");
    }
}

#[test]
fn a_refused_event_is_a_warning_and_the_file_copy_survives() {
    let dir = workdir("refused");
    let config = small_config(&dir, 20);
    let out = dir.join("out");
    let (url, _) = responder(2, 400);
    let o = cli(&config, &out, &[("OPENLINEAGE_URL", &url)]);
    assert!(o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("HTTP 400"));
    assert_eq!(events(&out).len(), 2);
}

#[test]
fn the_validator_rejects_broken_events() {
    let dir = workdir("negative");
    let config = small_config(&dir, 20);
    let out = dir.join("out");
    assert!(cli(&config, &out, &[]).status.success());
    let good = events(&out).pop().unwrap();
    let event_schema = "https://openlineage.io/spec/2-0-2/OpenLineage.json#/$defs/RunEvent";
    let mut no_run_id = good.clone();
    no_run_id["run"].as_object_mut().unwrap().remove("runId");
    assert!(
        !validate(event_schema, &no_run_id).is_empty(),
        "missing runId must fail"
    );
    let mut bad_type = good.clone();
    bad_type["eventType"] = json!("FINISHED");
    assert!(
        !validate(event_schema, &bad_type).is_empty(),
        "unknown eventType must fail"
    );
    let facet = &good["run"]["facets"]["stratigraphy_run"];
    let url = facet["_schemaURL"].as_str().unwrap();
    let mut bad_facet = facet.clone();
    bad_facet["config_source"] = json!("guess");
    assert!(
        !validate(url, &bad_facet).is_empty(),
        "custom facet enum must be enforced"
    );
    let mut no_producer = facet.clone();
    no_producer.as_object_mut().unwrap().remove("_producer");
    assert!(
        !validate(url, &no_producer).is_empty(),
        "custom facet must inherit RunFacet's required fields"
    );
}
