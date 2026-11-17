//! strat-forward — run the forward stratigraphic model from a JSON config.
//!
//!   strat-forward example-config > run.json
//!   strat-forward run --out DIR [--config run.json] [--deck] [--run-id UUID] [--openlineage-url URL]
//!
//! Without `--config` the built-in example runs. `--deck` also writes SECTION.DATA, an Eclipse-format
//! black-oil deck of the best sand window (default reservoir options), and deck_summary.json.
//!
//! Every run writes `prov.jsonld` (W3C PROV-O provenance of the stratigraphic record) and
//! `openlineage-events.json` (OpenLineage START and COMPLETE/FAIL events). With `--openlineage-url` or
//! `OPENLINEAGE_URL` set, the events are also POSTed there. `--run-id` sets the run's UUID (default: random);
//! it identifies the run in both documents.
//!
//! Outputs are staged and replace `DIR` whole only when the run succeeds (see `rundir`). A failed run leaves
//! `DIR` untouched and keeps what it wrote, with its FAIL event, in `DIR.failed-<run id>`. `DIR` must be
//! missing, empty, or a previous strat-forward output directory.

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use strat_forward::lineage::{self, FileRef, RunDescription, Sinks};
use strat_forward::prov::{self, RunInfo};
use strat_forward::rundir::RunDir;
use strat_forward::{export, Config, Model};
use serde_json::json;

const USAGE: &str = "usage: strat-forward run --out DIR [--config FILE] [--deck] [--run-id UUID] [--openlineage-url URL]\n       strat-forward example-config";

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

struct RunArgs {
    config: Option<PathBuf>,
    out: PathBuf,
    deck: bool,
    run_id: String,
    openlineage_url: Option<String>,
}

fn parse_run(args: &[String]) -> Result<RunArgs, String> {
    let (mut config, mut out, mut deck, mut run_id, mut url) = (None, None, false, None, None);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |flag: &str| {
            it.next()
                .cloned()
                .ok_or(format!("{flag} needs a value\n{USAGE}"))
        };
        match a.as_str() {
            "--config" => config = Some(PathBuf::from(value("--config")?)),
            "--out" => out = Some(PathBuf::from(value("--out")?)),
            "--deck" => deck = true,
            "--run-id" => run_id = Some(value("--run-id")?.to_ascii_lowercase()),
            "--openlineage-url" => url = Some(value("--openlineage-url")?),
            other => return Err(format!("unknown argument {other}\n{USAGE}")),
        }
    }
    let out = out.ok_or(format!("--out is required\n{USAGE}"))?;
    let run_id = run_id.unwrap_or_else(lineage::new_run_id);
    if !lineage::is_uuid(&run_id) {
        return Err(format!("--run-id must be a UUID, got {run_id}"));
    }
    Ok(RunArgs {
        config,
        out,
        deck,
        run_id,
        openlineage_url: url,
    })
}

fn real_main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("example-config") if args.len() == 1 => {
            let text =
                serde_json::to_string_pretty(&Config::example()).map_err(|e| e.to_string())?;
            println!("{text}");
            Ok(())
        }
        Some("run") => run(parse_run(&args[1..])?),
        _ => Err(USAGE.into()),
    }
}

fn run(a: RunArgs) -> Result<(), String> {
    let dir = RunDir::prepare(&a.out, &a.run_id).map_err(|e| e.to_string())?;
    let sinks = Sinks::new(
        dir.path().join("openlineage-events.json"),
        a.openlineage_url.clone(),
    );
    let started_at = lineage::now_rfc3339();
    let input = a
        .config
        .as_deref()
        .map(FileRef::of)
        .transpose()
        .map_err(|e| e.to_string())?;
    let source = if input.is_some() {
        "file"
    } else {
        "built_in_example"
    };
    let software = json!({ "name": env!("CARGO_PKG_NAME"), "version": env!("CARGO_PKG_VERSION") });
    let mut warnings: Vec<String> = Vec::new();
    // A lost HTTP delivery is a warning (the file copy exists); a lost file copy fails the run.
    let publish = |event: serde_json::Value, warnings: &mut Vec<String>| -> Result<(), String> {
        match sinks.publish(&event) {
            Ok(()) => Ok(()),
            Err(lineage::PublishError::Http(e)) => {
                warnings.push(e);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    };
    let mut desc = RunDescription {
        run_id: a.run_id.clone(),
        event_time: started_at.clone(),
        input: input.clone(),
        outputs: Vec::new(),
        facet: json!({ "software": software, "config_source": source }),
        error: None,
    };
    publish(lineage::run_event("START", &desc), &mut warnings)?;

    let result = forward(&a, dir.path(), &started_at, input.as_ref()).and_then(|(files, facet)| {
        // Outputs are hashed where they are, and reported where they will be after the commit.
        let outputs = files
            .iter()
            .map(|f| FileRef::of_as(f, &dir.target().join(f.file_name().unwrap_or_default())))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| e.to_string())?;
        Ok((outputs, facet))
    });
    desc.event_time = lineage::now_rfc3339();
    let outcome = match result {
        Ok((outputs, mut facet)) => {
            facet["software"] = software;
            facet["config_source"] = json!(source);
            desc.outputs = outputs;
            desc.facet = facet;
            publish(lineage::run_event("COMPLETE", &desc), &mut warnings)?;
            let names: Vec<String> = desc.outputs.iter().map(|o| o.path.clone()).collect();
            let done = dir.commit().map_err(|e| e.to_string())?;
            for n in names {
                println!("  wrote {n}");
            }
            println!(
                "  wrote {}",
                done.path.join("openlineage-events.json").display()
            );
            println!("  run id {}", a.run_id);
            if !done.atomic {
                warnings.push(format!(
                    "this filesystem has no atomic directory exchange; {} was replaced in two renames",
                    done.path.display()
                ));
            }
            if let Some(prev) = done.previous {
                println!("  previous outputs kept in {}", prev.display());
            }
            Ok(())
        }
        Err(e) => {
            desc.error = Some(e.clone());
            publish(lineage::run_event("FAIL", &desc), &mut warnings)?;
            let kept = dir.abandon().map_err(|err| {
                format!("{e}; and the failed run's files could not be kept: {err}")
            })?;
            Err(format!(
                "{e}\n  previous outputs untouched; this run's files and FAIL event are in {}",
                kept.display()
            ))
        }
    };
    for w in &warnings {
        eprintln!("warning: {w}");
    }
    outcome
}

/// The run itself, writing into `out`. Returns the files written and the body of the `stratigraphy_run` facet.
fn forward(
    a: &RunArgs,
    out: &Path,
    started_at: &str,
    input: Option<&FileRef>,
) -> Result<(Vec<PathBuf>, serde_json::Value), String> {
    let cfg: Config = match &a.config {
        Some(p) => {
            let raw =
                std::fs::read_to_string(p).map_err(|e| format!("read {}: {e}", p.display()))?;
            serde_json::from_str(&raw).map_err(|e| format!("config {}: {e}", p.display()))?
        }
        None => Config::example(),
    };
    let started = std::time::Instant::now();
    let mut model = Model::new(cfg).map_err(|e| e.to_string())?;
    let history = model.run().map_err(|e| e.to_string())?;
    let elapsed = started.elapsed();
    let mut files = export::write_all(out, &model, &history).map_err(|e| e.to_string())?;
    if a.deck {
        let d = strat_forward::reservoir::deck(&model, &Default::default())
            .map_err(|e| e.to_string())?;
        files.push(
            export::write_with(out, "SECTION.DATA", |w| w.write_all(d.text.as_bytes()))
                .map_err(|e| e.to_string())?,
        );
        files.push(
            export::write_json(out, "deck_summary.json", &d.summary, true)
                .map_err(|e| e.to_string())?,
        );
    }
    let info = RunInfo {
        run_id: a.run_id.clone(),
        started_at: started_at.to_string(),
        ended_at: lineage::now_rfc3339(),
        config: input.map(|f| (f.path.clone(), f.sha256.clone())),
    };
    files.push(
        export::write_with(out, "prov.jsonld", |w| prov::write_prov(&model, &info, w))
            .map_err(|e| e.to_string())?,
    );

    let budget = model.budget();
    let mass = export::max_relative_mass_error(&budget);
    let gap = export::surface_reconstruction_gap(&model, &history);
    let replay = strat_forward::record::replay(&model).map_err(|e| e.to_string())?;
    let cfg = model.cfg();
    println!(
        "ran {} record steps ({} yr substep) in {:.2?}",
        cfg.n_records(),
        model.substep_yr(),
        elapsed
    );
    for b in &budget {
        println!(
            "  {:<6} supplied {:>12.1} m²  stored {:>12.1} m²  outflow {:>10.1} m²  error {:+.3e} m²",
            b.class, b.supplied_m2, b.stored_m2, b.outflow_m2, b.error_m2
        );
    }
    println!("  max relative mass error {mass:.3e}");
    println!(
        "  surface-only reconstruction gap {gap:.3} m (record keeps what sampled surfaces miss)"
    );
    let facet = json!({
        "record_steps": cfg.n_records(),
        "cells": cfg.nx,
        "substep_yr": model.substep_yr(),
        "max_relative_mass_error": mass,
        "ledger_replay_max_diff_m": replay,
        "surface_reconstruction_gap_m": gap,
        "ledger_entries": model.ledger().iter().map(|o| o.len()).sum::<usize>(),
        "surviving_units": model.columns().iter().map(|c| c.len()).sum::<usize>(),
        "recycled_fraction": export::recycled_fraction(&model),
        "prov_document": "prov.jsonld",
    });
    Ok((files, facet))
}
