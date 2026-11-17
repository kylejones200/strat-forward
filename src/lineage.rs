//! OpenLineage RunEvents for a forward run.
//!
//! One run = one OpenLineage `Run`; `runId` is a UUID and is also the id of the run's PROV document
//! (`prov.jsonld`). The job is `strat-forward/strat-forward.run`. The input is the
//! config file (none for the built-in example); outputs are the files the run wrote. Events follow
//! OpenLineage spec 2-0-2 and the standard facets vendored in tests/fixtures/openlineage, which the tests
//! validate against, plus the custom facet in schemas/openlineage/StratigraphyFacets.json.
//!
//! Events always go to `openlineage-events.json` in the output directory. When `OPENLINEAGE_URL` is set
//! they are also sent with the standard OpenLineage HTTP transport (POST <url>/api/v1/lineage, with
//! `Authorization: Bearer $OPENLINEAGE_API_KEY` when that is set). No backend is named here.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde_json::{json, Map, Value};
use sha2::{Digest, Sha256};

pub const SCHEMA_URL: &str = "https://openlineage.io/spec/2-0-2/OpenLineage.json#/$defs/RunEvent";
pub const PRODUCER: &str = "https://github.com/kylejones200/strat-forward";
pub const JOB_NAMESPACE: &str = "strat-forward";
pub const JOB_NAME: &str = "strat-forward.run";
/// Schema of the custom facets (schemas/openlineage in this repo).
pub const FACET_SCHEMA: &str =
    "https://github.com/kylejones200/strat-forward/schemas/openlineage/StratigraphyFacets.json";

const PROCESSING_ENGINE: &str =
    "https://openlineage.io/spec/facets/1-1-1/ProcessingEngineRunFacet.json#/$defs/ProcessingEngineRunFacet";
const JOB_TYPE: &str =
    "https://openlineage.io/spec/facets/2-0-4/JobTypeJobFacet.json#/$defs/JobTypeJobFacet";
const DATA_SOURCE: &str =
    "https://openlineage.io/spec/facets/1-0-1/DatasourceDatasetFacet.json#/$defs/DatasourceDatasetFacet";
const DATASET_VERSION: &str =
    "https://openlineage.io/spec/facets/1-0-1/DatasetVersionDatasetFacet.json#/$defs/DatasetVersionDatasetFacet";
const OUTPUT_STATISTICS: &str = "https://openlineage.io/spec/facets/1-0-2/OutputStatisticsOutputDatasetFacet.json#/$defs/OutputStatisticsOutputDatasetFacet";
const ERROR_MESSAGE: &str =
    "https://openlineage.io/spec/facets/1-0-1/ErrorMessageRunFacet.json#/$defs/ErrorMessageRunFacet";

/// A file the run read or wrote.
#[derive(Debug, Clone)]
pub struct FileRef {
    pub path: String,
    pub sha256: String,
    pub bytes: u64,
}

impl FileRef {
    /// Hash and size of the file at `path`, reported under its own absolute path.
    pub fn of(path: &Path) -> crate::Result<Self> {
        Self::of_as(path, path)
    }

    /// Hash and size of the file at `on_disk`, reported as `reported` (for example its final location once
    /// a staged output directory is committed).
    pub fn of_as(on_disk: &Path, reported: &Path) -> crate::Result<Self> {
        let io = |e| crate::Error::io(on_disk, e);
        let mut f = std::fs::File::open(on_disk).map_err(io)?;
        let mut hasher = Sha256::new();
        let mut buf = vec![0u8; 1 << 16];
        let mut bytes = 0u64;
        loop {
            let n = f.read(&mut buf).map_err(io)?;
            if n == 0 {
                break;
            }
            hasher.update(&buf[..n]);
            bytes += n as u64;
        }
        let sha256 = hasher
            .finalize()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        Ok(FileRef {
            path: absolute(&reported.display().to_string()),
            sha256,
            bytes,
        })
    }
}

/// Everything a RunEvent says about the run.
#[derive(Debug, Clone)]
pub struct RunDescription {
    pub run_id: String,
    pub event_time: String,
    pub input: Option<FileRef>,
    pub outputs: Vec<FileRef>,
    /// Body of the `stratigraphy_run` facet.
    pub facet: Value,
    pub error: Option<String>,
}

/// An OpenLineage (namespace, name) for a local file: `file` + absolute path, as the naming spec says.
pub fn dataset_name(path: &str) -> (String, String) {
    ("file".to_string(), absolute(path))
}

/// A local path as a `file://` URI of its absolute path, percent-encoded outside the unreserved set.
pub fn location_uri(path: &str) -> String {
    let mut uri = String::from("file://");
    for b in absolute(path).bytes() {
        if b.is_ascii_alphanumeric() || b"-._~/".contains(&b) {
            uri.push(b as char);
        } else {
            uri.push_str(&format!("%{b:02X}"));
        }
    }
    uri
}

fn absolute(path: &str) -> String {
    let p = Path::new(path);
    if p.is_absolute() {
        return path.to_string();
    }
    std::env::current_dir()
        .map(|cwd| cwd.join(p).display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

fn facet(schema: &str, body: Value) -> Value {
    let mut obj = Map::new();
    obj.insert("_producer".into(), json!(PRODUCER));
    obj.insert("_schemaURL".into(), json!(schema));
    if let Value::Object(fields) = body {
        obj.extend(fields);
    }
    Value::Object(obj)
}

fn version_facet(sha: &str) -> Value {
    facet(
        DATASET_VERSION,
        json!({ "datasetVersion": format!("sha256:{sha}") }),
    )
}

fn input_dataset(f: &FileRef) -> Value {
    let (namespace, name) = dataset_name(&f.path);
    json!({
        "namespace": namespace,
        "name": name,
        "facets": {
            "dataSource": facet(DATA_SOURCE, json!({ "name": namespace, "uri": location_uri(&f.path) })),
            "version": version_facet(&f.sha256),
        }
    })
}

fn output_dataset(f: &FileRef) -> Value {
    let (namespace, name) = dataset_name(&f.path);
    json!({
        "namespace": namespace,
        "name": name,
        "facets": { "version": version_facet(&f.sha256) },
        "outputFacets": { "outputStatistics": facet(OUTPUT_STATISTICS, json!({ "size": f.bytes })) }
    })
}

/// `event_type` is START, COMPLETE or FAIL.
pub fn run_event(event_type: &str, run: &RunDescription) -> Value {
    let mut run_facets = Map::new();
    run_facets.insert(
        "processing_engine".into(),
        facet(
            PROCESSING_ENGINE,
            json!({ "name": env!("CARGO_PKG_NAME"), "version": env!("CARGO_PKG_VERSION") }),
        ),
    );
    run_facets.insert(
        "stratigraphy_run".into(),
        facet(
            &format!("{FACET_SCHEMA}#/$defs/StratigraphyRunFacet"),
            run.facet.clone(),
        ),
    );
    if let Some(message) = &run.error {
        run_facets.insert(
            "errorMessage".into(),
            facet(
                ERROR_MESSAGE,
                json!({ "message": message, "programmingLanguage": "rust" }),
            ),
        );
    }
    json!({
        "eventType": event_type,
        "eventTime": run.event_time,
        "producer": PRODUCER,
        "schemaURL": SCHEMA_URL,
        "run": { "runId": run.run_id, "facets": run_facets },
        "job": {
            "namespace": JOB_NAMESPACE,
            "name": JOB_NAME,
            "facets": {
                "jobType": facet(JOB_TYPE, json!({
                    "processingType": "BATCH",
                    "integration": "STRAT_FORWARD",
                    "jobType": "FORWARD_MODEL"
                }))
            }
        },
        "inputs": run.input.iter().map(input_dataset).collect::<Vec<_>>(),
        "outputs": run.outputs.iter().map(output_dataset).collect::<Vec<_>>(),
    })
}

/// A new random (v4) run id.
pub fn new_run_id() -> String {
    uuid::Uuid::new_v4().to_string()
}

/// True when `s` is a UUID in canonical form.
pub fn is_uuid(s: &str) -> bool {
    uuid::Uuid::parse_str(s).is_ok_and(|u| u.hyphenated().to_string() == s.to_ascii_lowercase())
}

/// Current UTC time as RFC 3339 with milliseconds.
pub fn now_rfc3339() -> String {
    let d = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = d.as_secs() as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // Civil date from days since 1970-01-01 (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + if month <= 2 { 1 } else { 0 };
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}.{:03}Z",
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60,
        d.subsec_millis()
    )
}

/// Where events go: always the evidence file, and the HTTP endpoint when one is configured.
pub struct Sinks {
    file: PathBuf,
    http: Option<HttpSink>,
}

impl Sinks {
    /// `url` overrides `OPENLINEAGE_URL`.
    pub fn new(file: PathBuf, url: Option<String>) -> Self {
        let url = url
            .or_else(|| std::env::var("OPENLINEAGE_URL").ok())
            .filter(|u| !u.trim().is_empty());
        let key = std::env::var("OPENLINEAGE_API_KEY")
            .ok()
            .filter(|k| !k.is_empty());
        Sinks {
            file,
            http: url.map(|u| HttpSink::new(&u, key)),
        }
    }

    /// Record the event: the file copy first, then the HTTP endpoint if any, so an HTTP failure never
    /// loses the evidence copy.
    pub fn publish(&self, event: &Value) -> Result<(), PublishError> {
        let file = |e: String| PublishError::File(format!("{}: {e}", self.file.display()));
        let mut events: Vec<Value> = match std::fs::read(&self.file) {
            Ok(bytes) => serde_json::from_slice(&bytes).map_err(|e| file(e.to_string()))?,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Vec::new(),
            Err(e) => return Err(file(e.to_string())),
        };
        events.push(event.clone());
        let body = serde_json::to_vec_pretty(&events).expect("serializable");
        std::fs::write(&self.file, body).map_err(|e| file(e.to_string()))?;
        match &self.http {
            Some(h) => h.publish(event).map_err(PublishError::Http),
            None => Ok(()),
        }
    }

    pub fn endpoint(&self) -> Option<&str> {
        self.http.as_ref().map(|h| h.endpoint.as_str())
    }
}

#[derive(Debug)]
pub enum PublishError {
    /// The evidence file could not be written.
    File(String),
    /// The file copy was written; the HTTP endpoint refused the event or could not be reached.
    Http(String),
}

impl std::fmt::Display for PublishError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PublishError::File(e) | PublishError::Http(e) => f.write_str(e),
        }
    }
}

/// OpenLineage HTTP transport: POST `<url>/api/v1/lineage`, bounded retries on network errors and 5xx.
struct HttpSink {
    endpoint: String,
    api_key: Option<String>,
    agent: ureq::Agent,
}

impl HttpSink {
    fn new(url: &str, api_key: Option<String>) -> Self {
        let base = url.trim().trim_end_matches('/');
        let endpoint = if base.ends_with("/api/v1/lineage") {
            base.to_string()
        } else {
            format!("{base}/api/v1/lineage")
        };
        let agent = ureq::AgentBuilder::new()
            .timeout_connect(Duration::from_secs(3))
            .timeout(Duration::from_secs(15))
            .build();
        HttpSink {
            endpoint,
            api_key,
            agent,
        }
    }

    fn publish(&self, event: &Value) -> Result<(), String> {
        let body = serde_json::to_string(event).expect("serializable");
        let mut last = String::new();
        for attempt in 0..3u32 {
            if attempt > 0 {
                std::thread::sleep(Duration::from_millis(200) * 2u32.pow(attempt - 1));
            }
            let mut req = self
                .agent
                .post(&self.endpoint)
                .set("content-type", "application/json");
            if let Some(key) = &self.api_key {
                req = req.set("authorization", &format!("Bearer {key}"));
            }
            match req.send_string(&body) {
                Ok(_) => return Ok(()),
                Err(ureq::Error::Status(status, resp)) if status < 500 => {
                    return Err(format!(
                        "OpenLineage event rejected by {} with HTTP {status}: {}",
                        self.endpoint,
                        resp.into_string().unwrap_or_default()
                    ))
                }
                Err(ureq::Error::Status(status, _)) => last = format!("HTTP {status}"),
                Err(ureq::Error::Transport(t)) => last = t.to_string(),
            }
        }
        Err(format!(
            "OpenLineage backend {} unavailable after 3 attempts: {last}",
            self.endpoint
        ))
    }
}
