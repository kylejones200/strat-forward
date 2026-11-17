//! W3C PROV-O provenance of the stratigraphic record, serialized as JSON-LD.
//!
//! The ledger maps onto PROV directly:
//! - every unit the run ever created is a `prov:Entity` (type `strat:Unit`);
//! - every ledger entry is a `prov:Activity` (type `strat:Deposit`, `strat:Erode` or `strat:Exchange`);
//! - a unit `prov:wasGeneratedBy` the entry that created it, and `prov:wasInfluencedBy` later entries that
//!   added to it;
//! - an entry that removed sediment from units `prov:used` them;
//! - a unit with nothing left at the end `prov:wasInvalidatedBy` the last entry that removed from it;
//! - the run is a `prov:Activity` `prov:wasAssociatedWith` the software agent; it `prov:used` the config.
//!
//! Model time is geological (Ma), which `xsd:dateTime` cannot express, so entries and units carry
//! `strat:startAgeMa` / `strat:endAgeMa` instead of `prov:startedAtTime`. Only the run itself has wall-clock
//! `prov:startedAtTime` / `prov:endedAtTime`.
//!
//! Surviving units also carry the origin tracer: the recycled fraction of their solids (sediment that had
//! been deposited and eroded before) and the mean age at which their sediment entered the basin.
//!
//! The `strat:` vocabulary is local to this crate (`urn:strat-forward:vocab#`); README.md defines
//! its terms. Node identifiers are URNs scoped by the run id, which is also the OpenLineage `runId`.

use std::collections::HashMap;

use serde_json::{json, Map, Value};

use crate::model::{Contact, Layer, Model, OpKind};

pub const PROV_NS: &str = "http://www.w3.org/ns/prov#";
pub const STRAT_NS: &str = "urn:strat-forward:vocab#";

/// What the run knows about itself, shared by the PROV document and the OpenLineage events.
#[derive(Clone, Debug)]
pub struct RunInfo {
    /// UUID; also the OpenLineage runId.
    pub run_id: String,
    /// RFC 3339 wall-clock times.
    pub started_at: String,
    pub ended_at: String,
    /// The config file read (path) with its SHA-256; None for the built-in example.
    pub config: Option<(String, String)>,
}

fn id(v: &str) -> Value {
    json!({ "@id": v })
}

fn ids(v: Vec<String>) -> Value {
    Value::Array(v.into_iter().map(|s| id(&s)).collect())
}

/// A property name for a per-class quantity: `strat:<prefix>_<class>_<suffix>`.
fn class_key(prefix: &str, class: &str, suffix: &str) -> String {
    let clean: String = class
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    format!("strat:{prefix}_{clean}_{suffix}")
}

struct UnitFacts {
    cell: usize,
    step: usize,
    t_first: f64,
    t_last: f64,
    deposited: Vec<f64>,
    generated_by: Option<String>,
    influenced_by: Vec<String>,
    last_removal: Option<String>,
    contact: Option<Contact>,
    initial: bool,
}

fn context() -> Value {
    json!({ "prov": PROV_NS, "xsd": "http://www.w3.org/2001/XMLSchema#", "strat": STRAT_NS })
}

/// The whole document as one JSON value. Holds every node in memory: use `write_prov` for real runs.
pub fn prov_document(model: &Model, run: &RunInfo) -> Value {
    let mut graph = Vec::new();
    generate(model, run, &mut |n| {
        graph.push(n);
        Ok(())
    })
    .expect("collecting into memory cannot fail");
    json!({ "@context": context(), "@graph": graph })
}

/// Stream the document as JSON-LD, one node at a time; memory does not grow with the number of nodes
/// beyond the per-unit facts the relations need.
pub fn write_prov(
    model: &Model,
    run: &RunInfo,
    w: &mut impl std::io::Write,
) -> std::io::Result<()> {
    write!(w, "{{\"@context\":")?;
    serde_json::to_writer(&mut *w, &context()).map_err(std::io::Error::other)?;
    write!(w, ",\"@graph\":[")?;
    let mut first = true;
    generate(model, run, &mut |n| {
        if !first {
            w.write_all(b",")?;
        }
        first = false;
        serde_json::to_writer(&mut *w, &n).map_err(std::io::Error::other)
    })?;
    w.write_all(b"]}")
}

fn generate(
    model: &Model,
    run: &RunInfo,
    emit: &mut dyn FnMut(Value) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let cfg = &model.cfg;
    let nc = cfg.classes.len();
    let base = format!("urn:strat-forward:{}", run.run_id);
    let run_iri = format!("urn:uuid:{}", run.run_id);
    let agent_iri = format!(
        "urn:strat-forward:software:{}",
        env!("CARGO_PKG_VERSION")
    );
    let unit_iri = |u: u64| format!("{base}:unit:{u}");
    let op_iri = |x: usize, k: usize| format!("{base}:op:{x}:{k}");
    let age = |t_yr: f64| cfg.start_age_ma - t_yr / 1e6;

    let surviving: HashMap<u64, &Layer> =
        model.columns.iter().flatten().map(|l| (l.id, l)).collect();
    let mut units: Vec<(u64, UnitFacts)> = Vec::new();
    let mut index: HashMap<u64, usize> = HashMap::new();
    for (x, init) in model.initial_units.iter().enumerate() {
        for (u, solid) in init {
            index.insert(*u, units.len());
            units.push((
                *u,
                UnitFacts {
                    cell: x,
                    step: 0,
                    t_first: 0.0,
                    t_last: 0.0,
                    deposited: solid.clone(),
                    generated_by: None,
                    influenced_by: Vec::new(),
                    last_removal: None,
                    contact: Some(Contact::Basement),
                    initial: true,
                },
            ));
        }
    }

    for (x, ops) in model.ledger.iter().enumerate() {
        for (k, op) in ops.iter().enumerate() {
            let iri = op_iri(x, k);
            if let Some(u) = op.into_unit {
                if op.created.is_some() {
                    index.insert(u, units.len());
                    units.push((
                        u,
                        UnitFacts {
                            cell: x,
                            step: op.step,
                            t_first: op.t_start_yr,
                            t_last: op.t_end_yr,
                            deposited: vec![0.0; nc],
                            generated_by: Some(iri.clone()),
                            influenced_by: Vec::new(),
                            last_removal: None,
                            contact: op.created.clone(),
                            initial: false,
                        },
                    ));
                }
                if let Some(&i) = index.get(&u) {
                    let f = &mut units[i].1;
                    for c in 0..nc {
                        f.deposited[c] += op.deposited[c];
                    }
                    f.t_last = f.t_last.max(op.t_end_yr);
                    if op.created.is_none() {
                        f.influenced_by.push(iri.clone());
                    }
                }
            }
            for r in &op.removed {
                if let Some(&i) = index.get(&r.unit) {
                    units[i].1.last_removal = Some(iri.clone());
                }
            }

            let kind = match op.kind {
                OpKind::Deposit => "strat:Deposit",
                OpKind::Erode => "strat:Erode",
                OpKind::Exchange => "strat:Exchange",
            };
            let mut node = Map::new();
            node.insert("@id".into(), json!(iri));
            node.insert("@type".into(), json!(["prov:Activity", kind]));
            node.insert("strat:run".into(), id(&run_iri));
            node.insert("strat:cell".into(), json!(x));
            node.insert("strat:recordStep".into(), json!(op.step));
            node.insert("strat:startAgeMa".into(), json!(age(op.t_start_yr)));
            node.insert("strat:endAgeMa".into(), json!(age(op.t_end_yr)));
            node.insert("strat:substeps".into(), json!(op.substeps));
            node.insert("strat:bypassM".into(), json!(op.bypass_m));
            for (c, cls) in cfg.classes.iter().enumerate() {
                if op.deposited[c] > 0.0 {
                    node.insert(
                        class_key("deposited", &cls.name, "m"),
                        json!(op.deposited[c]),
                    );
                }
                let removed: f64 = op.removed.iter().map(|r| r.solid[c]).sum();
                if removed > 0.0 {
                    node.insert(class_key("removed", &cls.name, "m"), json!(removed));
                }
            }
            if !op.removed.is_empty() {
                let mut units_used: Vec<u64> = op.removed.iter().map(|r| r.unit).collect();
                units_used.sort_unstable();
                units_used.dedup();
                let used: Vec<String> = units_used.into_iter().map(unit_iri).collect();
                node.insert("prov:used".into(), ids(used));
            }
            emit(Value::Object(node))?;
        }
    }

    for (u, f) in &units {
        let mut node = Map::new();
        node.insert("@id".into(), json!(unit_iri(*u)));
        node.insert("@type".into(), json!(["prov:Entity", "strat:Unit"]));
        node.insert("strat:run".into(), id(&run_iri));
        node.insert("strat:cell".into(), json!(f.cell));
        node.insert(
            "strat:xCenterM".into(),
            json!((f.cell as f64 + 0.5) * model.dx),
        );
        node.insert("strat:recordStep".into(), json!(f.step));
        node.insert("strat:startAgeMa".into(), json!(age(f.t_first)));
        node.insert("strat:endAgeMa".into(), json!(age(f.t_last)));
        if f.initial {
            node.insert("strat:initialSediment".into(), json!(true));
        }
        for (c, cls) in cfg.classes.iter().enumerate() {
            node.insert(
                class_key("deposited", &cls.name, "m"),
                json!(f.deposited[c]),
            );
        }
        match &f.contact {
            Some(Contact::Basement) => {
                node.insert("strat:baseContact".into(), json!("basement"));
            }
            Some(Contact::Conformable) => {
                node.insert("strat:baseContact".into(), json!("conformable"));
            }
            Some(Contact::Diastem { missing_steps }) => {
                node.insert("strat:baseContact".into(), json!("diastem"));
                node.insert("strat:missingSteps".into(), json!(missing_steps));
            }
            Some(Contact::Erosional {
                cut_m,
                truncated_units,
            }) => {
                node.insert("strat:baseContact".into(), json!("erosional"));
                node.insert("strat:cutM".into(), json!(cut_m));
                node.insert(
                    "strat:truncates".into(),
                    ids(truncated_units.iter().map(|t| unit_iri(*t)).collect()),
                );
            }
            None => {}
        }
        if let Some(g) = &f.generated_by {
            node.insert("prov:wasGeneratedBy".into(), id(g));
        }
        if !f.influenced_by.is_empty() {
            node.insert("prov:wasInfluencedBy".into(), ids(f.influenced_by.clone()));
        }
        match surviving.get(u) {
            Some(l) => {
                for (c, cls) in cfg.classes.iter().enumerate() {
                    node.insert(class_key("surviving", &cls.name, "m"), json!(l.solid[c]));
                }
                node.insert(
                    "strat:recycledFraction".into(),
                    json!(l.recycled_fraction()),
                );
                if l.total_solid() > 0.0 {
                    let entry = l.mean_entry_yr();
                    node.insert("strat:meanEntryAgeMa".into(), json!(age(entry)));
                    // Mean time from entering the basin to the middle of the unit's deposition.
                    let mid = 0.5 * (l.t_first_yr + l.t_last_yr);
                    node.insert(
                        "strat:meanTransitMyr".into(),
                        json!(((mid - entry) / 1e6).max(0.0)),
                    );
                }
            }
            None => {
                if let Some(r) = &f.last_removal {
                    node.insert("prov:wasInvalidatedBy".into(), id(r));
                }
            }
        }
        emit(Value::Object(node))?;
    }

    let mut run_node = json!({
        "@id": run_iri,
        "@type": ["prov:Activity", "strat:ForwardRun"],
        "prov:startedAtTime": { "@value": run.started_at, "@type": "xsd:dateTime" },
        "prov:endedAtTime": { "@value": run.ended_at, "@type": "xsd:dateTime" },
        "prov:wasAssociatedWith": id(&agent_iri),
        "strat:startAgeMa": cfg.start_age_ma,
        "strat:endAgeMa": age(cfg.duration_yr),
        "strat:recordSteps": cfg.n_records(),
        "strat:cells": cfg.nx,
    });
    if let Some((path, sha)) = &run.config {
        let cfg_iri = format!("{base}:config");
        run_node["prov:used"] = id(&cfg_iri);
        emit(json!({
            "@id": cfg_iri,
            "@type": ["prov:Entity", "strat:Config"],
            "prov:atLocation": id(&crate::lineage::location_uri(path)),
            "strat:sha256": sha,
        }))?;
    }
    emit(run_node)?;
    emit(json!({
        "@id": agent_iri,
        "@type": ["prov:Agent", "prov:SoftwareAgent"],
        "strat:name": env!("CARGO_PKG_NAME"),
        "strat:version": env!("CARGO_PKG_VERSION"),
    }))
}
