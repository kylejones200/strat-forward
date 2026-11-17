//! Sensitivity and convergence study on the reference configuration.
//!
//!   cargo run --release --example sensitivity [-- OUT_DIR]
//!
//! Families (each compared with its finest member):
//! - time step: cfl 0.4, 0.2, 0.1
//! - layer clock: layer_dt 20, 10, 5, 2 kyr (record interval fixed at 20 kyr)
//! - grid: nx 100, 200, 400
//! - shoreline width: 0, 2, 5, 10 m, each at nx 100, 200, 400
//! - active layer: 0.5, 1, 2 m
//!
//! Writes `sensitivity.md` and `sensitivity.json`.

use std::sync::Mutex;
use std::time::Instant;

use strat_forward::model::Contact;
use strat_forward::{export, record, Config, Model};
use serde::Serialize;

#[derive(Clone, Serialize)]
struct Row {
    family: String,
    label: String,
    reference: String,
    nx: usize,
    cfl: f64,
    layer_dt_yr: f64,
    shoreline_width_m: f64,
    active_layer_m: f64,
    substep_yr: f64,
    runtime_s: f64,
    final_shoreline_km: f64,
    max_shoreline_km: f64,
    gross_eroded_m2: f64,
    gross_deposited_m2: f64,
    erosional_contacts: usize,
    max_cut_m: f64,
    proximal_sand_fraction: f64,
    surface_gap_m: f64,
    mass_error: f64,
    /// RMS difference of preserved horizons (every 100 kyr, sampled every 1 km) against the reference run.
    horizon_rms_vs_ref_m: Option<f64>,
    #[serde(skip)]
    horizons: Vec<f64>,
}

struct Job {
    family: &'static str,
    label: String,
    reference: String,
    cfg: Config,
}

fn sample_horizons(m: &Model, s: &record::Stratigraphy) -> Vec<f64> {
    // Horizons at record steps 5, 10, … (every 100 kyr at 20 kyr records), linear in x at 1 km spacing.
    let nx = m.cfg().nx;
    let nt = m.cfg().n_records();
    let mut out = Vec::new();
    for k in (5..=nt).step_by(5) {
        for km in 1..(m.cfg().length_m / 1000.0) as usize {
            let x = km as f64 * 1000.0 / m.dx() - 0.5;
            let j0 = (x.floor().max(0.0) as usize).min(nx - 2);
            let f = (x - j0 as f64).clamp(0.0, 1.0);
            let a = s.preserved_horizon_m[k * nx + j0];
            let b = s.preserved_horizon_m[k * nx + j0 + 1];
            out.push(a + (b - a) * f);
        }
    }
    out
}

fn run(job: &Job) -> Row {
    let started = Instant::now();
    let mut m = Model::new(job.cfg.clone()).expect("config");
    let h = m.run().unwrap();
    let runtime_s = started.elapsed().as_secs_f64();
    let s = record::stratigraphy(&m);
    let dx = m.dx();
    let (dep, _) = record::gross(&m);
    let mut max_cut = 0.0f64;
    for col in &s.columns {
        for u in col {
            if let Contact::Erosional { cut_m, .. } = u.base {
                max_cut = max_cut.max(cut_m);
            }
        }
    }
    let nx = m.cfg().nx;
    let (mut sand, mut tot) = (0.0, 0.0);
    for col in m.columns().iter().take(nx / 4) {
        for l in col {
            sand += l.solid[0];
            tot += l.total_solid();
        }
    }
    let shore: Vec<f64> = h
        .shoreline_x_m
        .iter()
        .map(|v| v.unwrap_or(m.cfg().length_m) / 1000.0)
        .collect();
    Row {
        family: job.family.to_string(),
        label: job.label.clone(),
        reference: job.reference.clone(),
        nx,
        cfl: m.cfg().transport.cfl,
        layer_dt_yr: m.cfg().layer_dt(),
        shoreline_width_m: m.cfg().transport.shoreline_width_m,
        active_layer_m: m.cfg().transport.active_layer_m,
        substep_yr: m.substep_yr(),
        runtime_s,
        final_shoreline_km: *shore.last().unwrap(),
        max_shoreline_km: shore.iter().cloned().fold(0.0, f64::max),
        gross_eroded_m2: m.eroded_m().iter().sum::<f64>() * dx,
        gross_deposited_m2: dep.iter().sum::<f64>() * dx,
        erosional_contacts: s.contacts.erosional,
        max_cut_m: max_cut,
        proximal_sand_fraction: sand / tot.max(1e-30),
        surface_gap_m: export::surface_reconstruction_gap(&m, &h),
        mass_error: export::max_relative_mass_error(&m.budget()),
        horizon_rms_vs_ref_m: None,
        horizons: sample_horizons(&m, &s),
    }
}

fn main() {
    let out_dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "out/sensitivity".into());
    let base = Config::example();
    let mut jobs: Vec<Job> = Vec::new();
    let mut add =
        |family: &'static str, label: String, reference: String, f: &dyn Fn(&mut Config)| {
            let mut c = base.clone();
            f(&mut c);
            jobs.push(Job {
                family,
                label,
                reference,
                cfg: c,
            });
        };
    for cfl in [0.4, 0.2, 0.1] {
        add("time step", format!("cfl {cfl}"), "cfl 0.1".into(), &|c| {
            c.transport.cfl = cfl
        });
    }
    for l in [20_000.0, 10_000.0, 5_000.0, 2_000.0] {
        add(
            "layer clock",
            format!("layer {} kyr", l / 1000.0),
            "layer 2 kyr".into(),
            &|c| c.layer_dt_yr = Some(l),
        );
    }
    for nx in [100, 200, 400] {
        add("grid", format!("nx {nx}"), "nx 400".into(), &|c| c.nx = nx);
    }
    for w in [0.0, 2.0, 5.0, 10.0] {
        for nx in [100, 200, 400] {
            add(
                "shoreline width",
                format!("w {w} m, nx {nx}"),
                format!("w {w} m, nx 400"),
                &|c| {
                    c.transport.shoreline_width_m = w;
                    c.nx = nx;
                },
            );
        }
    }
    for a in [0.5, 1.0, 2.0] {
        add(
            "active layer",
            format!("active {a} m"),
            "active 0.5 m".into(),
            &|c| c.transport.active_layer_m = a,
        );
    }

    // Optional second argument: comma-separated family names to run, e.g. "layer clock,active layer".
    if let Some(only) = std::env::args().nth(2) {
        let keep: Vec<String> = only.split(',').map(|f| f.trim().to_string()).collect();
        jobs.retain(|j| keep.iter().any(|k| k == j.family));
    }
    // Longest first, on all cores.
    jobs.sort_by_key(|j| {
        std::cmp::Reverse(j.cfg.nx * 1000 / (j.cfg.transport.cfl * 100.0) as usize)
    });
    let n = jobs.len();
    let queue = Mutex::new(jobs.into_iter().enumerate().collect::<Vec<_>>());
    let results = Mutex::new(Vec::new());
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    eprintln!("{n} runs on {threads} threads");
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let next = queue.lock().unwrap().pop();
                let Some((i, job)) = next else { break };
                let row = run(&job);
                eprintln!(
                    "  done {:<16} {:<22} {:>6.1}s",
                    row.family, row.label, row.runtime_s
                );
                results.lock().unwrap().push((i, row));
            });
        }
    });
    let mut rows: Vec<Row> = results
        .into_inner()
        .unwrap()
        .into_iter()
        .map(|(_, r)| r)
        .collect();
    let order = [
        "time step",
        "layer clock",
        "grid",
        "shoreline width",
        "active layer",
    ];
    rows.sort_by(|a, b| {
        let fa = order.iter().position(|f| *f == a.family).unwrap();
        let fb = order.iter().position(|f| *f == b.family).unwrap();
        (fa, &a.label).cmp(&(fb, &b.label))
    });
    let snapshot = rows.clone();
    for r in rows.iter_mut() {
        if let Some(reference) = snapshot
            .iter()
            .find(|x| x.family == r.family && x.label == r.reference)
        {
            let n = r.horizons.len().min(reference.horizons.len());
            let ss: f64 = (0..n)
                .map(|k| (r.horizons[k] - reference.horizons[k]).powi(2))
                .sum();
            r.horizon_rms_vs_ref_m = Some((ss / n as f64).sqrt());
        }
    }

    let mut md = String::from(
        "# Sensitivity study\n\nDemo case (`Config::example()`), one parameter family at a time. Horizon RMS compares preserved horizons every 100 kyr, sampled every 1 km, against the reference run named in the table.\n\n",
    );
    md.push_str(&format!(
        "Base config: max_unit_solid_m = {:?}, layer_dt_yr = {:?}, record_dt_yr = {}.\n\n",
        base.max_unit_solid_m, base.layer_dt_yr, base.record_dt_yr
    ));
    md.push_str("| family | run | substep yr | runtime s | final shoreline km | max shoreline km | gross eroded m² | gross deposited m² | erosional contacts | max cut m | proximal sand | surface gap m | horizon RMS vs ref m | mass error |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    for r in &rows {
        md.push_str(&format!(
            "| {} | {} | {:.2} | {:.1} | {:.2} | {:.2} | {:.0} | {:.0} | {} | {:.2} | {:.3} | {:.3} | {} | {:.0e} |\n",
            r.family,
            r.label,
            r.substep_yr,
            r.runtime_s,
            r.final_shoreline_km,
            r.max_shoreline_km,
            r.gross_eroded_m2,
            r.gross_deposited_m2,
            r.erosional_contacts,
            r.max_cut_m,
            r.proximal_sand_fraction,
            r.surface_gap_m,
            r.horizon_rms_vs_ref_m.map_or("–".into(), |v| format!("{v:.3}")),
            r.mass_error
        ));
    }
    std::fs::create_dir_all(&out_dir).expect("out dir");
    std::fs::write(format!("{out_dir}/sensitivity.md"), &md).expect("write md");
    std::fs::write(
        format!("{out_dir}/sensitivity.json"),
        serde_json::to_string_pretty(&rows).unwrap(),
    )
    .expect("write json");
    println!("{md}");
}
