//! Origin tracer: where each unit's sediment came from (fresh supply or recycled from eroded units) and when
//! it entered the basin. The tracer is passive, so these tests check its own bookkeeping: exact closure,
//! physical bounds, and that it responds to erosion.

use strat_forward::config::{Basement, Boundary};
use strat_forward::{Config, Model};

/// The benchmark ramp with one sea-level cycle: erosion during the fall, an open distal end.
fn cycle() -> Config {
    let mut c = Config::example();
    c.length_m = 40_000.0;
    c.nx = 80;
    c.duration_yr = 400_000.0;
    c.record_dt_yr = 10_000.0;
    c.initial.basement = Basement::Ramp {
        top_m: 10.0,
        bottom_m: -200.0,
    };
    c.supply.flux_m2_per_yr = 1.0;
    c.subsidence.proximal_m_per_yr = 0.0;
    c.subsidence.distal_m_per_yr = 0.0;
    c.sea_level.amplitude_m = 20.0;
    c.sea_level.period_yr = 400_000.0;
    c.sea_level.trend_m_per_yr = 0.0;
    c.transport.distal_boundary = Boundary::Open;
    c
}

fn run(c: Config) -> Model {
    let mut m = Model::new(c).unwrap();
    m.run().unwrap();
    m
}

fn stored(
    m: &Model,
    f: impl Fn(&strat_forward::model::Layer, usize) -> f64,
    i: usize,
) -> f64 {
    m.columns().iter().flatten().map(|l| f(l, i)).sum::<f64>() * m.dx()
}

#[test]
fn the_tracer_budget_closes_exactly() {
    for c in [cycle(), {
        let mut c = cycle();
        c.initial.sediment_solid_m = 2.0;
        c
    }] {
        let m = run(c);
        let b = &m.tracer();
        for i in 0..m.cfg().classes.len() {
            let entry = stored(&m, |l, i| l.entry_moment[i], i) + b.outflow_entry_m2yr[i];
            let rel = (entry - b.supplied_entry_m2yr[i]).abs() / b.supplied_entry_m2yr[i];
            assert!(rel < 1e-10, "class {i}: entry moment closes to {rel:e}");
            let rec = stored(&m, |l, i| l.recycled[i], i) + b.outflow_recycled_m2[i];
            let rel = (rec - b.recycled_created_m2[i]).abs() / b.recycled_created_m2[i].max(1e-300);
            assert!(rel < 1e-9, "class {i}: recycled closes to {rel:e}");
            assert!(b.recycled_created_m2[i] > 0.0, "the cycle erodes");
        }
    }
}

#[test]
fn every_unit_is_a_mix_of_its_sources() {
    let m = run(cycle());
    let tol = 1e-9;
    for col in m.columns() {
        for l in col {
            for i in 0..l.solid.len() {
                assert!(
                    l.recycled[i] >= -tol && l.recycled[i] <= l.solid[i] + tol,
                    "unit {} class {i}",
                    l.id
                );
                // Sediment cannot have entered after the unit finished forming, nor before the run began.
                let s = l.solid[i];
                if s > 1e-9 {
                    let mean = l.entry_moment[i] / s;
                    assert!(
                        mean >= -1e-6 && mean <= l.t_last_yr + 1e-6,
                        "unit {} class {i}: entry {mean}",
                        l.id
                    );
                }
            }
        }
    }
}

#[test]
fn erosion_is_what_makes_sediment_recycled() {
    let mut steady = cycle();
    steady.sea_level.amplitude_m = 0.0;
    let rec = |m: &Model| {
        m.columns()
            .iter()
            .flatten()
            .map(|l| l.recycled.iter().sum::<f64>())
            .sum::<f64>()
    };
    let tot = |m: &Model| {
        m.columns()
            .iter()
            .flatten()
            .map(|l| l.total_solid())
            .sum::<f64>()
    };
    let (a, b) = (run(steady), run(cycle()));
    let (fa, fb) = (rec(&a) / tot(&a), rec(&b) / tot(&b));
    println!("recycled fraction: constant sea level {fa:.4}, cycle {fb:.4}");
    assert!(
        fb > 2.0 * fa,
        "a sea-level fall recycles more ({fb} vs {fa})"
    );
}

// ---- W3C PROV-O document ----

use strat_forward::prov::{prov_document, RunInfo};
use serde_json::Value;
use std::collections::{HashMap, HashSet};

fn doc(m: &Model) -> Value {
    let info = RunInfo {
        run_id: "00000000-0000-4000-8000-000000000001".into(),
        started_at: "2026-11-09T00:00:00.000Z".into(),
        ended_at: "2026-11-09T00:01:00.000Z".into(),
        config: None,
    };
    prov_document(m, &info)
}

fn types(n: &Value) -> Vec<&str> {
    n["@type"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t.as_str().unwrap())
        .collect()
}

/// Every `{"@id": …}` reference inside a node, with the property it sits under.
fn refs(n: &Value) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for (k, v) in n.as_object().unwrap() {
        let items: Vec<&Value> = match v {
            Value::Array(a) => a.iter().collect(),
            other => vec![other],
        };
        for it in items {
            if let Some(id) = it.get("@id").and_then(Value::as_str) {
                if k != "@id" && !id.starts_with("file://") {
                    out.push((k.clone(), id.to_string()));
                }
            }
        }
    }
    out
}

#[test]
fn the_prov_document_is_closed_typed_and_uses_the_w3c_namespace() {
    let m = run(cycle());
    let d = doc(&m);
    assert_eq!(d["@context"]["prov"], "http://www.w3.org/ns/prov#");
    let graph = d["@graph"].as_array().unwrap();
    let mut defined = HashSet::new();
    for n in graph {
        let id = n["@id"].as_str().unwrap();
        assert!(defined.insert(id.to_string()), "{id} defined twice");
        let t = types(n);
        // Entity, Activity and Agent are disjoint here.
        let kinds = ["prov:Entity", "prov:Activity", "prov:Agent"]
            .iter()
            .filter(|k| t.contains(k))
            .count();
        assert_eq!(kinds, 1, "{id} has types {t:?}");
    }
    for n in graph {
        for (k, r) in refs(n) {
            assert!(defined.contains(&r), "{} {k} -> undefined {r}", n["@id"]);
        }
    }
    // One activity per ledger entry plus the run; one entity per unit ever created.
    let activities = graph
        .iter()
        .filter(|n| types(n).contains(&"prov:Activity"))
        .count();
    assert_eq!(
        activities,
        m.ledger().iter().map(|o| o.len()).sum::<usize>() + 1
    );
    let created = m
        .ledger()
        .iter()
        .flatten()
        .filter(|o| o.created.is_some())
        .count();
    let initial: usize = m.initial_units().iter().map(|u| u.len()).sum();
    let units = graph
        .iter()
        .filter(|n| types(n).contains(&"strat:Unit"))
        .count();
    assert_eq!(units, created + initial);
}

#[test]
fn prov_generation_use_and_invalidation_follow_the_ledger_in_time() {
    let m = run(cycle());
    let d = doc(&m);
    let graph = d["@graph"].as_array().unwrap();
    let node: HashMap<&str, &Value> = graph
        .iter()
        .map(|n| (n["@id"].as_str().unwrap(), n))
        .collect();
    let start_age = |id: &str| node[id]["strat:startAgeMa"].as_f64().unwrap();
    let surviving: HashSet<u64> = m.columns().iter().flatten().map(|l| l.id).collect();
    // Units used by each activity.
    let mut used_by: HashMap<String, Vec<String>> = HashMap::new();
    for n in graph.iter().filter(|n| types(n).contains(&"prov:Activity")) {
        for (k, r) in refs(n) {
            if k == "prov:used" && r.contains(":unit:") {
                used_by
                    .entry(r)
                    .or_default()
                    .push(n["@id"].as_str().unwrap().to_string());
            }
        }
    }
    let (mut invalidated, mut generated) = (0, 0);
    for n in graph.iter().filter(|n| types(n).contains(&"strat:Unit")) {
        let id = n["@id"].as_str().unwrap();
        let unit: u64 = id.rsplit(':').next().unwrap().parse().unwrap();
        let gen = n.get("prov:wasGeneratedBy").and_then(|g| g["@id"].as_str());
        let inv = n
            .get("prov:wasInvalidatedBy")
            .and_then(|g| g["@id"].as_str());
        // PROV-CONSTRAINTS: at most one generation and one invalidation (single-valued here by construction);
        // surviving units are never invalidated, and every unit that is gone was invalidated.
        assert_eq!(inv.is_some(), !surviving.contains(&unit), "{id}");
        // Ages decrease forward in time: generation precedes use, use precedes invalidation.
        if let Some(g) = gen {
            generated += 1;
            for u in used_by.get(id).into_iter().flatten() {
                assert!(start_age(g) >= start_age(u), "{id}: used before generated");
            }
            if let Some(i) = inv {
                assert!(
                    start_age(g) >= start_age(i),
                    "{id}: invalidated before generated"
                );
            }
        }
        if let Some(i) = inv {
            invalidated += 1;
            for u in used_by.get(id).into_iter().flatten() {
                assert!(
                    start_age(u) >= start_age(i),
                    "{id}: used after invalidation"
                );
            }
            assert!(
                used_by.get(id).is_some_and(|v| v.iter().any(|u| u == i)),
                "{id}: invalidated by a non-user"
            );
        }
    }
    println!("units generated {generated}, invalidated {invalidated}");
    assert!(invalidated > 0, "the sea-level fall removes whole units");
}

#[test]
fn surviving_units_carry_their_origin_in_prov() {
    let m = run(cycle());
    let d = doc(&m);
    let start = m.cfg().start_age_ma;
    for n in d["@graph"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|n| types(n).contains(&"strat:Unit"))
    {
        if n.get("prov:wasInvalidatedBy").is_some() {
            continue;
        }
        let r = n["strat:recycledFraction"].as_f64().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&r));
        if let Some(entry) = n.get("strat:meanEntryAgeMa").and_then(Value::as_f64) {
            // Sediment entered no earlier than the run start and no later than the unit's last deposit.
            assert!(entry <= start + 1e-9 && entry >= n["strat:endAgeMa"].as_f64().unwrap() - 1e-9);
        }
    }
}
