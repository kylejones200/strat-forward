//! Numerical acceptance tests. These check the solver, not the geology.

use strat_forward::config::{Basement, Boundary, GrainClass};
use strat_forward::export::{max_relative_mass_error, surface_reconstruction_gap};
use strat_forward::record::{self, ChronoState};
use strat_forward::{Config, Model};

fn one_class(phi0: f64) -> Vec<GrainClass> {
    vec![GrainClass {
        name: "sed".into(),
        land_k_mult: 1.0,
        marine_k_mult: 1.0,
        phi0,
        compaction_length_m: 2000.0,
    }]
}

fn small(cfg: &mut Config) {
    cfg.length_m = 20_000.0;
    cfg.nx = 40;
    cfg.duration_yr = 200_000.0;
    cfg.record_dt_yr = 10_000.0;
    cfg.transport.k_land_m2_per_yr = 2_000.0;
    cfg.transport.k_marine_m2_per_yr = 400.0;
}

#[test]
fn closed_basin_stores_exactly_what_was_supplied() {
    let mut cfg = Config::example();
    small(&mut cfg);
    cfg.transport.distal_boundary = Boundary::Closed;
    let mut m = Model::new(cfg.clone()).unwrap();
    m.run().unwrap();
    let b = m.budget();
    for c in &b {
        assert_eq!(c.outflow_m2, 0.0);
    }
    let supplied: f64 = b.iter().map(|c| c.supplied_m2).sum();
    let stored: f64 = b.iter().map(|c| c.stored_m2).sum();
    assert!((supplied - cfg.supply.flux_m2_per_yr * cfg.duration_yr).abs() < 1e-6 * supplied);
    assert!(
        (stored - supplied).abs() < 1e-9 * supplied,
        "stored {stored} supplied {supplied}"
    );
    assert!(max_relative_mass_error(&b) < 1e-10);
}

#[test]
fn open_basin_balances_storage_and_outflow() {
    let mut cfg = Config::example();
    small(&mut cfg);
    cfg.transport.distal_boundary = Boundary::Open;
    cfg.transport.k_deep_m2_per_yr = 200.0; // push material to the distal edge
    let mut m = Model::new(cfg).unwrap();
    m.run().unwrap();
    let b = m.budget();
    let out: f64 = b.iter().map(|c| c.outflow_m2).sum();
    assert!(out > 0.0, "expected some outflow");
    assert!(max_relative_mass_error(&b) < 1e-10, "{b:?}");
}

#[test]
fn zero_flux_boundaries_conserve_initial_sediment() {
    let mut cfg = Config::example();
    small(&mut cfg);
    cfg.supply.flux_m2_per_yr = 0.0;
    cfg.initial.sediment_solid_m = 30.0;
    cfg.transport.distal_boundary = Boundary::Closed;
    let mut m = Model::new(cfg.clone()).unwrap();
    m.run().unwrap();
    let b = m.budget();
    let stored: f64 = b.iter().map(|c| c.stored_m2).sum();
    let initial = 30.0 * cfg.length_m;
    assert!(
        (stored - initial).abs() < 1e-9 * initial,
        "stored {stored} initial {initial}"
    );
    // And the sediment actually moved: the ramp diffused.
    let first = m.columns()[0].iter().map(|l| l.total_solid()).sum::<f64>();
    assert!((first - 30.0).abs() > 1e-3);
}

#[test]
fn prescribed_subsidence_moves_the_basement_exactly() {
    let mut cfg = Config::example();
    small(&mut cfg);
    cfg.supply.flux_m2_per_yr = 0.0;
    let mut m = Model::new(cfg.clone()).unwrap();
    let b0 = m.basement().to_vec();
    m.run().unwrap();
    for (j, b) in b0.iter().enumerate() {
        let rate = cfg.subsidence.rate((j as f64 + 0.5) / cfg.nx as f64);
        let expect = b - rate * cfg.duration_yr;
        // ~16k substeps of accumulated subtraction: allow rounding relative to the elevation.
        let err = (m.basement()[j] - expect).abs();
        assert!(err < 1e-12 * expect.abs() * 1e3, "cell {j}: err {err}");
        assert!(
            (m.eta()[j] - m.basement()[j]).abs() < 1e-12,
            "bare basement is the surface"
        );
    }
}

#[test]
fn pure_diffusion_matches_the_analytic_decay_of_a_cosine_mode() {
    let (l, nx, k, amp) = (10_000.0_f64, 100usize, 100.0_f64, 2.0_f64);
    let dx = l / nx as f64;
    let mut cfg = Config::example();
    cfg.length_m = l;
    cfg.nx = nx;
    cfg.duration_yr = 100_000.0;
    cfg.record_dt_yr = 10_000.0;
    cfg.classes = one_class(0.0);
    cfg.supply.flux_m2_per_yr = 0.0;
    cfg.supply.fractions = vec![1.0];
    cfg.subsidence.proximal_m_per_yr = 0.0;
    cfg.subsidence.distal_m_per_yr = 0.0;
    cfg.sea_level.base_m = -1e6; // everything is land: uniform K
    cfg.sea_level.amplitude_m = 0.0;
    cfg.sea_level.trend_m_per_yr = 0.0;
    cfg.transport.k_land_m2_per_yr = k;
    cfg.transport.distal_boundary = Boundary::Closed;
    cfg.initial.sediment_solid_m = 50.0;
    cfg.initial.basement = Basement::Profile(
        (0..nx)
            .map(|j| amp * (std::f64::consts::PI * (j as f64 + 0.5) * dx / l).cos())
            .collect(),
    );
    let mut m = Model::new(cfg.clone()).unwrap();
    let a0 = (m.eta()[0] - m.eta()[nx - 1]) / 2.0;
    m.run().unwrap();
    let a1 = (m.eta()[0] - m.eta()[nx - 1]) / 2.0;
    let ratio = a1 / a0;
    let continuous = (-k * (std::f64::consts::PI / l).powi(2) * cfg.duration_yr).exp();
    let discrete_lambda = 2.0 / (dx * dx) * (1.0 - (std::f64::consts::PI * dx / l).cos());
    let discrete = (-k * discrete_lambda * cfg.duration_yr).exp();
    assert!(
        (ratio - discrete).abs() < 2e-4,
        "ratio {ratio} discrete {discrete}"
    );
    assert!(
        (ratio - continuous).abs() < 2e-3,
        "ratio {ratio} continuous {continuous}"
    );
    // Shape stays a cosine: midpoint stays at the mean.
    let mean: f64 = m.eta().iter().sum::<f64>() / nx as f64;
    assert!((m.eta()[nx / 2] + m.eta()[nx / 2 - 1] - 2.0 * mean).abs() < 1e-6);
}

fn gap_for(cfg: Config) -> f64 {
    let mut m = Model::new(cfg).unwrap();
    let h = m.run().unwrap();
    surface_reconstruction_gap(&m, &h)
}

#[test]
fn the_surface_reconstruction_gap_never_exceeds_what_was_eroded() {
    // Surfaces sampled at record times miss cut-and-fill inside a step. They cannot miss more than was cut.
    let mut cfg = Config::example();
    small(&mut cfg);
    let mut m = Model::new(cfg).unwrap();
    let h = m.run().unwrap();
    let gap = surface_reconstruction_gap(&m, &h);
    let cut = m.eroded_m().iter().cloned().fold(0.0, f64::max);
    assert!(gap <= cut + 1e-9, "gap {gap} cut {cut}");
}

#[test]
fn with_no_erosion_the_surface_grid_suffix_minimum_is_exactly_the_record() {
    // One class, one diffusivity everywhere (all marine, no depth decay), steady supply into a closed basin.
    let mut cfg = Config::example();
    small(&mut cfg);
    cfg.classes = one_class(0.5);
    cfg.supply.fractions = vec![1.0];
    cfg.sea_level.base_m = 1e6;
    cfg.sea_level.amplitude_m = 0.0;
    cfg.sea_level.trend_m_per_yr = 0.0;
    cfg.transport.depth_decay_m = 1e12;
    let mut m = Model::new(cfg).unwrap();
    let h = m.run().unwrap();
    let cut = m.eroded_m().iter().cloned().fold(0.0, f64::max);
    let gap = surface_reconstruction_gap(&m, &h);
    assert!(cut < 1e-9, "setup eroded {cut}");
    assert!(gap < 1e-9, "gap {gap}");
}

#[test]
fn the_surface_reconstruction_gap_comes_from_sampling_and_shrinks_with_it() {
    // A surface sampled every record step cannot see cut-and-fill inside the step; the event record can.
    let mut coarse = Config::example();
    small(&mut coarse);
    coarse.classes = one_class(0.5);
    coarse.supply.fractions = vec![1.0];
    let mut fine = coarse.clone();
    fine.record_dt_yr = coarse.record_dt_yr / 4.0;
    let (g_coarse, g_fine) = (gap_for(coarse), gap_for(fine));
    assert!(g_coarse > 0.0);
    assert!(g_fine < g_coarse, "fine {g_fine} coarse {g_coarse}");
}

#[test]
fn occupancy_uses_the_four_decision_states() {
    let tol = 0.5;
    assert_eq!(record::classify(4.0, 4.0, tol), ChronoState::Deposition);
    assert_eq!(record::classify(0.2, 0.0, tol), ChronoState::Hiatus);
    assert_eq!(record::classify(-3.0, 0.0, tol), ChronoState::Erosion);
    assert_eq!(record::classify(6.0, 0.0, tol), ChronoState::Vacuity);
    assert_eq!(ChronoState::Vacuity.decision(), "rise_removed");
    assert_eq!(ChronoState::Deposition.code(), 1);
}

#[test]
fn the_example_progrades_and_partitions_sand_landward() {
    let mut m = Model::new(Config::example()).unwrap();
    let h = m.run().unwrap();
    let first = h.shoreline_x_m.iter().flatten().next().copied().unwrap();
    let last = h.shoreline_x_m.last().copied().flatten().unwrap();
    assert!(last > first + 5_000.0, "shoreline {first} -> {last}");
    let sand_frac = |range: std::ops::Range<usize>| {
        let (mut s, mut t) = (0.0, 0.0);
        for j in range {
            for l in &m.columns()[j] {
                s += l.solid[0];
                t += l.total_solid();
            }
        }
        s / t
    };
    let nx = m.cfg().nx;
    assert!(
        sand_frac(0..nx / 4) > sand_frac(3 * nx / 4..nx),
        "proximal should be sandier"
    );
}

#[test]
fn a_smooth_shoreline_removes_spurious_erosion_under_steady_rise() {
    // A sharp land/marine diffusivity switch makes shoreline cells cut and refill even while sea level
    // only rises. Blending it over a few meters of water depth cuts that erosion by ~6 orders of magnitude.
    // The ~1e-5 m left is time-discretization error in the distal cells, not the shoreline: a 4x smaller time
    // step cuts it 16-fold (9.9e-6 to 6.2e-7 m) and each grid doubling halves it.
    let run = |width: f64| {
        let mut cfg = Config::example();
        small(&mut cfg);
        cfg.duration_yr = 30_000.0;
        cfg.record_dt_yr = 100.0;
        cfg.sea_level.amplitude_m = 0.0;
        cfg.classes = one_class(0.5);
        cfg.supply.fractions = vec![1.0];
        cfg.transport.shoreline_width_m = width;
        let mut m = Model::new(cfg).unwrap();
        let h = m.run().unwrap();
        let cut = m.eroded_m().iter().cloned().fold(0.0, f64::max);
        (cut, surface_reconstruction_gap(&m, &h))
    };
    let (sharp_cut, _) = run(0.0);
    let (smooth_cut, smooth_gap) = run(5.0);
    assert!(
        sharp_cut > 1.0,
        "sharp switch should show the artifact, cut {sharp_cut}"
    );
    assert!(smooth_cut < 1e-3, "smooth blend still erodes {smooth_cut}");
    assert!(smooth_gap < 1e-3, "gap {smooth_gap}");
}

#[test]
fn configs_without_a_shoreline_width_keep_the_sharp_switch() {
    let mut v = serde_json::to_value(Config::example()).unwrap();
    v["transport"]
        .as_object_mut()
        .unwrap()
        .remove("shoreline_width_m");
    let cfg: Config = serde_json::from_value(v).unwrap();
    assert_eq!(cfg.transport.shoreline_width_m, 0.0);
}

#[test]
fn rounding_level_changes_do_not_get_ledger_entries_of_their_own() {
    // In bypass cells supply and outflow cancel to ~1e-17 m per substep. Such changes join the previous
    // entry of the same step; as entries of their own they flipped kind every substep (24% of the reference ledger).
    let mut c: Config =
        serde_json::from_str(&std::fs::read_to_string("examples/prograding_wedge.json").unwrap())
            .unwrap();
    c.duration_yr = 200_000.0;
    let mut m = Model::new(c).unwrap();
    m.run().unwrap();
    let moved = |o: &strat_forward::model::Op| {
        let dep: f64 = o.deposited.iter().sum();
        let rem: f64 = o.removed.iter().flat_map(|r| r.solid.iter()).sum();
        dep.max(rem)
    };
    let mut stray = 0;
    for ops in m.ledger() {
        for w in ops.windows(2) {
            if w[1].step == w[0].step && w[1].created.is_none() && moved(&w[1]) < 1e-12 {
                stray += 1;
            }
        }
    }
    assert_eq!(
        stray, 0,
        "rounding-level entries after the first of their step"
    );
    assert!(record::replay(&m).unwrap() < 1e-9);
}

#[test]
fn the_unit_cap_is_on_unless_a_config_sets_it_to_null() {
    use strat_forward::config::DEFAULT_MAX_UNIT_SOLID_M;
    let mut v = serde_json::to_value(Config::example()).unwrap();
    let obj = v.as_object_mut().unwrap();
    // Missing: the default cap.
    obj.remove("max_unit_solid_m");
    let missing: Config = serde_json::from_value(v.clone()).unwrap();
    assert_eq!(missing.max_unit_solid_m, Some(DEFAULT_MAX_UNIT_SOLID_M));
    // null: no cap.
    v["max_unit_solid_m"] = serde_json::Value::Null;
    let off: Config = serde_json::from_value(v.clone()).unwrap();
    assert_eq!(off.max_unit_solid_m, None);
    // A value: that value.
    v["max_unit_solid_m"] = serde_json::json!(0.3);
    let set: Config = serde_json::from_value(v).unwrap();
    assert_eq!(set.max_unit_solid_m, Some(0.3));
    // Written out and read back, each keeps its meaning: off is written as null, not left out.
    for c in [missing, off, set] {
        let text = serde_json::to_string(&c).unwrap();
        assert!(text.contains("\"max_unit_solid_m\":"), "{text}");
        let back: Config = serde_json::from_str(&text).unwrap();
        assert_eq!(back.max_unit_solid_m, c.max_unit_solid_m);
    }
}
