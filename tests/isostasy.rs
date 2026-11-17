//! Flexural isostasy: the plate solver against exact results, then the coupled model.

use strat_forward::isostasy::{deflection, flexural_parameter, G};

const RHO_M: f64 = 3300.0;
const RHO_W: f64 = 1025.0;

/// Rigidity for elastic thickness `te` with E = 70 GPa, ν = 0.25.
fn rigidity(te: f64) -> f64 {
    7.0e10 * te.powi(3) / (12.0 * (1.0 - 0.0625))
}

#[test]
fn a_uniform_load_gives_the_airy_deflection_whatever_the_rigidity() {
    // Tolerance 1e-7: a stiff plate on a fine grid is ill-conditioned (bending / restoring ~ 1e7 per cell),
    // so double precision gives ~1e-9 relative here, not 1e-15.
    let n = 101;
    let k = vec![(RHO_M - RHO_W) * G; n];
    let q = vec![2.0e6; n];
    for te in [0.0, 5e3, 30e3, 80e3] {
        let w = deflection(1000.0, rigidity(te), &k, &q);
        let airy = 2.0e6 / ((RHO_M - RHO_W) * G);
        for (j, v) in w.iter().enumerate() {
            assert!(
                (v - airy).abs() < 1e-7 * airy,
                "te {te} cell {j}: {v} vs {airy}"
            );
        }
    }
}

#[test]
fn zero_rigidity_is_local_isostasy_cell_by_cell() {
    let n = 50;
    let k: Vec<f64> = (0..n)
        .map(|j| {
            if j < 20 {
                RHO_M * G
            } else {
                (RHO_M - RHO_W) * G
            }
        })
        .collect();
    let q: Vec<f64> = (0..n).map(|j| 1e5 * (j as f64).sin()).collect();
    let w = deflection(500.0, 0.0, &k, &q);
    for j in 0..n {
        assert!((w[j] - q[j] / k[j]).abs() < 1e-12);
    }
}

#[test]
fn total_force_and_moment_balance_exactly() {
    // Arbitrary load and restoring force: with free ends the plate carries no net force or moment.
    // Exact in exact arithmetic; 1e-7 for the same conditioning reason as the uniform-load test.
    let n = 120;
    let dx = 750.0;
    let k: Vec<f64> = (0..n)
        .map(|j| {
            if j % 7 < 3 {
                RHO_M * G
            } else {
                (RHO_M - RHO_W) * G
            }
        })
        .collect();
    let q: Vec<f64> = (0..n)
        .map(|j| 3e6 * (-((j as f64 - 30.0) / 8.0).powi(2)).exp())
        .collect();
    let w = deflection(dx, rigidity(25e3), &k, &q);
    let force: f64 = (0..n).map(|j| k[j] * w[j] - q[j]).sum();
    let moment: f64 = (0..n).map(|j| j as f64 * (k[j] * w[j] - q[j])).sum();
    let scale: f64 = q.iter().sum();
    assert!(force.abs() < 1e-7 * scale, "force residual {force}");
    assert!(
        moment.abs() < 1e-7 * scale * n as f64,
        "moment residual {moment}"
    );
}

/// Line load V0 (N/m) on an infinite plate (Turcotte & Schubert eq. 3-127):
/// w(x) = w0 e^(−x/α) (cos x/α + sin x/α), w0 = V0 α³ / (8 D).
fn line_load_error(dx: f64) -> (f64, f64, f64, f64) {
    let d = rigidity(20e3);
    let restoring = (RHO_M - RHO_W) * G;
    let alpha = flexural_parameter(d, restoring);
    // Plate 40 α long, so the free ends are far from the load.
    let n = (40.0 * alpha / dx) as usize | 1;
    let mid = n / 2;
    let v0 = 1.0e12;
    let mut q = vec![0.0; n];
    q[mid] = v0 / dx;
    let w = deflection(dx, d, &vec![restoring; n], &q);
    let w0 = v0 * alpha.powi(3) / (8.0 * d);
    let exact = |x: f64| {
        let s = x.abs() / alpha;
        w0 * (-s).exp() * (s.cos() + s.sin())
    };
    let err = (0..n)
        .map(|j| (w[j] - exact((j as f64 - mid as f64) * dx)).abs())
        .fold(0.0, f64::max);
    // Forebulge: the most negative deflection, expected at x = π α with height −w0 e^(−π).
    let (jb, wb) = (mid..n)
        .map(|j| (j, w[j]))
        .fold((mid, 0.0), |a, b| if b.1 < a.1 { b } else { a });
    let xb = (jb - mid) as f64 * dx;
    (
        err / w0,
        (xb - std::f64::consts::PI * alpha).abs() / alpha,
        wb / w0,
        alpha,
    )
}

#[test]
fn a_line_load_matches_the_analytic_flexure_profile_and_converges() {
    let (e1, xb1, wb1, alpha) = line_load_error(2000.0);
    let (e2, xb2, wb2, _) = line_load_error(1000.0);
    let (e3, _, _, _) = line_load_error(500.0);
    println!("alpha {alpha:.0} m; max error / w0: {e1:.2e} {e2:.2e} {e3:.2e}");
    assert!(e2 < 0.01, "1 km grid: {e2}");
    // Second order: halving dx cuts the error about fourfold.
    assert!(e1 / e2 > 3.0 && e2 / e3 > 3.0, "{e1} {e2} {e3}");
    // Forebulge at π α (within one cell) with height −e^(−π) w0 = −0.0432 w0.
    for (xb, wb) in [(xb1, wb1), (xb2, wb2)] {
        assert!(xb < 2000.0 / alpha + 1e-9, "forebulge offset {xb} α");
        assert!(
            (wb + (-std::f64::consts::PI).exp()).abs() < 2e-3,
            "forebulge {wb}"
        );
    }
}

// ---- Coupled model ----

use strat_forward::config::{Basement, Boundary, Isostasy};
use strat_forward::{record, Config, Model};

const RHO_G: f64 = 2650.0;

/// A deep basin that stays under water: constant sea level, basement 150–400 m deep, modest supply.
fn deep(te: f64, response_time_yr: Option<f64>) -> Config {
    let mut c = Config::example();
    c.length_m = 60_000.0;
    c.nx = 60;
    c.duration_yr = 200_000.0;
    c.record_dt_yr = 20_000.0;
    c.initial.basement = Basement::Ramp {
        top_m: -150.0,
        bottom_m: -400.0,
    };
    c.sea_level.amplitude_m = 0.0;
    c.sea_level.trend_m_per_yr = 0.0;
    c.supply.flux_m2_per_yr = 2.0;
    c.transport.distal_boundary = Boundary::Open;
    c.isostasy = Some(Isostasy {
        elastic_thickness_m: te,
        mantle_density_kg_m3: RHO_M,
        water_density_kg_m3: RHO_W,
        grain_density_kg_m3: RHO_G,
        youngs_modulus_pa: 7.0e10,
        poisson_ratio: 0.25,
        response_time_yr,
    });
    c
}

fn run(c: Config) -> Model {
    let mut m = Model::new(c).unwrap();
    m.run().unwrap();
    m
}

/// Solid added to each column since the start.
fn added_solid(m: &Model) -> Vec<f64> {
    let s0 = m.cfg().initial.sediment_solid_m;
    m.columns()
        .iter()
        .map(|c| c.iter().map(|l| l.total_solid()).sum::<f64>() - s0)
        .collect()
}

#[test]
fn local_isostasy_in_the_model_is_exact_while_the_basin_stays_under_water() {
    let m = run(deep(0.0, None));
    let sea = m.cfg().sea_level.at(m.t_yr());
    assert!(
        m.eta().iter().all(|&z| z < sea),
        "the basin must stay submerged"
    );
    let w = m.isostatic_deflection().unwrap();
    let ratio = (RHO_G - RHO_W) / (RHO_M - RHO_W);
    for (j, s) in added_solid(&m).iter().enumerate() {
        let expect = s * ratio;
        assert!(
            (w[j] - expect).abs() < 1e-9 * expect.abs().max(1.0),
            "cell {j}: {} vs {expect}",
            w[j]
        );
    }
    assert!(
        w.iter().cloned().fold(0.0, f64::max) > 1.0,
        "the test must load the plate"
    );
}

#[test]
fn a_flexed_basin_balances_total_force_and_moment_against_its_sediment_load() {
    // For any rigidity, total buoyancy equals total excess sediment weight, and so do their moments.
    let m = run(deep(20e3, None));
    let w = m.isostatic_deflection().unwrap();
    let s = added_solid(&m);
    let force_w: f64 = w.iter().sum::<f64>() * (RHO_M - RHO_W);
    let force_s: f64 = s.iter().sum::<f64>() * (RHO_G - RHO_W);
    let moment_w: f64 =
        w.iter().enumerate().map(|(j, v)| j as f64 * v).sum::<f64>() * (RHO_M - RHO_W);
    let moment_s: f64 =
        s.iter().enumerate().map(|(j, v)| j as f64 * v).sum::<f64>() * (RHO_G - RHO_W);
    assert!(
        (force_w - force_s).abs() < 1e-7 * force_s,
        "{force_w} vs {force_s}"
    );
    assert!(
        (moment_w - moment_s).abs() < 1e-7 * moment_s,
        "{moment_w} vs {moment_s}"
    );
    // A stiff plate spreads the load: the deflection is smoother than the local one.
    let local: Vec<f64> = s
        .iter()
        .map(|v| v * (RHO_G - RHO_W) / (RHO_M - RHO_W))
        .collect();
    let rough = |v: &[f64]| {
        v.windows(3)
            .map(|x| (x[0] - 2.0 * x[1] + x[2]).abs())
            .sum::<f64>()
    };
    assert!(
        rough(w) < rough(&local),
        "{} vs {}",
        rough(w),
        rough(&local)
    );
}

#[test]
fn isostasy_keeps_conservation_and_the_ledger_exact() {
    for te in [0.0, 20e3] {
        let mut c = deep(te, None);
        // Include land and a sea-level cycle, so cells cross the shoreline.
        c.initial.basement = Basement::Ramp {
            top_m: 15.0,
            bottom_m: -250.0,
        };
        c.sea_level.amplitude_m = 20.0;
        c.sea_level.period_yr = 100_000.0;
        let m = run(c);
        let mass = strat_forward::export::max_relative_mass_error(&m.budget());
        assert!(mass < 1e-9, "te {te}: mass error {mass}");
        let replay = record::replay(&m).unwrap();
        assert!(replay < 1e-6, "te {te}: replay {replay}");
    }
}

#[test]
fn the_response_time_delays_the_deflection() {
    let eq = run(deep(20e3, None));
    let slow = run(deep(20e3, Some(1.0e9)));
    let fast = run(deep(20e3, Some(1.0)));
    let max = |m: &Model| {
        m.isostatic_deflection()
            .unwrap()
            .iter()
            .cloned()
            .fold(0.0, f64::max)
    };
    assert!(
        max(&slow) < 1e-3 * max(&eq),
        "slow {} vs eq {}",
        max(&slow),
        max(&eq)
    );
    assert!(
        (max(&fast) - max(&eq)).abs() < 1e-3 * max(&eq),
        "fast {} vs eq {}",
        max(&fast),
        max(&eq)
    );
}

#[test]
fn isostasy_off_leaves_the_config_and_basement_untouched() {
    let mut c = deep(0.0, None);
    c.isostasy = None;
    assert!(!serde_json::to_string(&c).unwrap().contains("isostasy"));
    let m = run(c);
    assert!(m.isostatic_deflection().is_none());
}

/// Consecutive single-substep ledger entries of different kinds: a cell flip-flopping every substep.
fn substep_flips(m: &Model) -> usize {
    m.ledger()
        .iter()
        .map(|ops| {
            ops.windows(2)
                .filter(|w| w[0].substeps == 1 && w[1].substeps == 1 && w[0].kind != w[1].kind)
                .count()
        })
        .sum()
}

#[test]
fn the_plate_equilibrium_is_continuous_at_the_shoreline() {
    // With tectonic subsidence, a cell that crosses sea level must see the same equilibrium from both sides.
    // A jump there (water in the tectonic subsidence counted on one side only) pushes shoreline cells back
    // across every substep: in the reference configuration it multiplied the ledger 20-fold.
    let base = |iso: bool| {
        let mut c: Config = serde_json::from_str(
            &std::fs::read_to_string("examples/prograding_wedge.json").unwrap(),
        )
        .unwrap();
        c.duration_yr = 400_000.0;
        if !iso {
            c.isostasy = None;
        } else {
            c.isostasy = deep(0.0, None).isostasy;
        }
        c
    };
    let off = substep_flips(&run(base(false)));
    let local = substep_flips(&run(base(true)));
    assert!(
        local < 2 * off,
        "local isostasy {local} flips vs {off} without"
    );
}
