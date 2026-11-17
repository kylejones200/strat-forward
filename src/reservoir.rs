//! Reservoir deck from the preserved stratigraphy, in ECLIPSE format (FIELD units), for
//! black-oil simulators that read Eclipse decks, such as OPM Flow.
//!
//! - Grid: corner point, NX × 1 × NZ over a window of the section. One layer per record step, youngest at
//!   the top, bounded by the preserved horizons. Layers thinner than `min_thickness_m` are inactive and
//!   `PINCH` connects across them.
//! - The preserved column is re-compacted with `overburden_m` of extra burial (default: the datum depth), so
//!   thickness and porosity match the depth the deck places them at. PORO: mean porosity of the step's
//!   surviving units; NTG: sand fraction of the solids.
//! - PERM: per class, Kozeny–Carman scaling of a reference permeability with porosity,
//!   k = k_ref · (φ/φ_ref)³ · ((1−φ_ref)/(1−φ))²; classes are mixed log-linearly by sand fraction.
//!   PERMZ = kv/kh · PERMX. These are a stand-in rock model, not calibrated.
//! - Fluids, rock and relative permeability: the SPE1 benchmark tables. Initialized overpressured
//!   (4800 psi at the top, as SPE1) so the oil is undersaturated; oil–water contact below the window,
//!   gas–oil contact above it.
//! - Wells: gas injector in the proximal column, oil producer in the distal column, one completion per
//!   active layer, SPE1-style controls.

use std::fmt::Write as _;

use serde::Serialize;

use crate::model::Model;

const FT: f64 = 3.280_839_895;

#[derive(Clone, Debug, Serialize)]
pub struct ReservoirOptions {
    /// Window along the section. None picks the `window_m`-wide window with the most sand pore volume.
    pub x_start_m: Option<f64>,
    pub window_m: f64,
    /// Depth of elevation 0 below the reservoir datum: depth = datum_depth_m − elevation.
    pub datum_depth_m: f64,
    /// Extra burial used to re-compact the column to reservoir depth. None means `datum_depth_m`.
    pub overburden_m: Option<f64>,
    /// Strike width of the 2D grid.
    pub width_m: f64,
    pub min_thickness_m: f64,
    pub sand_k_ref_md: f64,
    pub mud_k_ref_md: f64,
    pub phi_ref: f64,
    pub kv_kh: f64,
    /// Report steps of 30 days.
    pub months: usize,
    pub producer_oil_rate_stb_d: f64,
    pub injector_gas_rate_mscf_d: f64,
}

impl Default for ReservoirOptions {
    fn default() -> Self {
        ReservoirOptions {
            x_start_m: None,
            window_m: 10_000.0,
            datum_depth_m: 2_500.0,
            overburden_m: None,
            width_m: 300.0,
            min_thickness_m: 0.1,
            sand_k_ref_md: 1_000.0,
            mud_k_ref_md: 0.01,
            phi_ref: 0.3,
            kv_kh: 0.1,
            months: 24,
            producer_oil_rate_stb_d: 2_000.0,
            injector_gas_rate_mscf_d: 10_000.0,
        }
    }
}

/// One cell of the stratigraphic grid, before it is written.
#[derive(Clone, Debug, Serialize)]
pub struct GridSummary {
    pub x_start_m: f64,
    pub x_end_m: f64,
    pub nx: usize,
    pub nz: usize,
    pub active_cells: usize,
    pub top_depth_ft: f64,
    pub bottom_depth_ft: f64,
    pub pore_volume_rb: f64,
    pub mean_poro_active: f64,
    pub mean_ntg_active: f64,
}

pub struct Deck {
    pub text: String,
    pub summary: GridSummary,
}

fn kozeny(k_ref: f64, phi: f64, phi_ref: f64) -> f64 {
    let phi = phi.clamp(1e-3, 0.95);
    k_ref * (phi / phi_ref).powi(3) * ((1.0 - phi_ref) / (1.0 - phi)).powi(2)
}

/// The preserved column re-compacted as if buried a further `overburden_m`: per (column, step) thickness,
/// porosity and sand fraction, and the horizon elevations those thicknesses imply.
struct Reburied {
    nx: usize,
    /// (thickness m, porosity, sand fraction), index `x * (nt + 1) + step`.
    props: Vec<(f64, f64, f64)>,
    /// Elevation of the top of all units with step <= t, index `t * nx + x`.
    horizon: Vec<f64>,
    nt: usize,
}

impl Reburied {
    fn new(model: &Model, overburden_m: f64) -> Self {
        let nx = model.cfg.nx;
        let nt = model.cfg.n_records();
        let classes = &model.cfg.classes;
        let mut props = vec![(0.0, 0.0, 0.0); nx * (nt + 1)];
        let mut horizon = vec![0.0; (nt + 1) * nx];
        for (x, col) in model.columns.iter().enumerate() {
            // Thickness of each unit at its burial depth plus the overburden, from the top down.
            let mut bulk = vec![0.0; col.len()];
            let mut z = overburden_m;
            for (k, l) in col.iter().enumerate().rev() {
                let tot = l.total_solid();
                if tot <= 0.0 {
                    continue;
                }
                let (mut phi0, mut c) = (0.0, 0.0);
                for (s, cls) in l.solid.iter().zip(classes) {
                    phi0 += s / tot * cls.phi0;
                    c += s / tot * cls.compaction_length_m;
                }
                bulk[k] = crate::compaction::bulk_thickness(tot, z, phi0, c);
                z += bulk[k];
            }
            let mut acc = vec![(0.0, 0.0, 0.0); nt + 1]; // bulk, solid, sand per step
            for (l, b) in col.iter().zip(&bulk) {
                let s = l.step.min(nt);
                acc[s].0 += b;
                acc[s].1 += l.total_solid();
                acc[s].2 += l.solid[0];
            }
            let mut top = model.basement[x];
            for t in 0..=nt {
                let (b, s, sand) = acc[t];
                top += b;
                horizon[t * nx + x] = top;
                props[x * (nt + 1) + t] = if b > 0.0 && s > 0.0 {
                    (b, 1.0 - s / b, sand / s)
                } else {
                    (0.0, 0.0, 0.0)
                };
            }
        }
        Reburied {
            nx,
            props,
            horizon,
            nt,
        }
    }
    fn at(&self, x: usize, step: usize) -> (f64, f64, f64) {
        self.props[x * (self.nt + 1) + step]
    }
    fn horizon(&self, x: usize, t: usize) -> f64 {
        self.horizon[t * self.nx + x]
    }
}

fn pick_window(model: &Model, rb: &Reburied, cells: usize) -> usize {
    let nx = model.cfg.nx;
    let nt = model.cfg.n_records();
    let score: Vec<f64> = (0..nx)
        .map(|x| {
            (1..=nt)
                .map(|t| {
                    let (h, phi, f) = rb.at(x, t);
                    h * phi * f
                })
                .sum()
        })
        .collect();
    let cells = cells.min(nx);
    (0..=nx - cells)
        .max_by(|&a, &b| {
            let s = |i: usize| score[i..i + cells].iter().sum::<f64>();
            s(a).total_cmp(&s(b))
        })
        .unwrap_or(0)
}

fn write_array(out: &mut String, name: &str, values: &[f64], digits: usize) {
    let _ = writeln!(out, "{name}");
    for chunk in values.chunks(8) {
        let line: Vec<String> = chunk.iter().map(|v| format!("{v:.digits$}")).collect();
        let _ = writeln!(out, "  {}", line.join(" "));
    }
    let _ = writeln!(out, "/");
}

pub fn deck(model: &Model, opt: &ReservoirOptions) -> crate::Result<Deck> {
    deck_inner(model, opt).map_err(crate::Error::Unsupported)
}

fn deck_inner(model: &Model, opt: &ReservoirOptions) -> Result<Deck, String> {
    let nt = model.cfg.n_records();
    let dx = model.dx;
    let rb = Reburied::new(model, opt.overburden_m.unwrap_or(opt.datum_depth_m));
    let cells = ((opt.window_m / dx).round() as usize).clamp(2, model.cfg.nx);
    let x0 = match opt.x_start_m {
        Some(x) => ((x / dx).floor() as usize).min(model.cfg.nx - cells),
        None => pick_window(model, &rb, cells),
    };
    let ni = cells;
    let nk = nt;
    // (preserved horizons come from the reburied column below)
    // Horizon elevation at column center: top of all units with step <= t.
    let horizon = |x: usize, t: usize| rb.horizon(x, t);
    let depth_ft = |elev: f64| (opt.datum_depth_m - elev) * FT;
    // Pillar (edge) depths: average of neighboring column centers.
    let edge = |e: usize, t: usize| -> f64 {
        let l = if e == 0 { x0 } else { x0 + e - 1 };
        let r = if e == ni { x0 + ni - 1 } else { x0 + e };
        depth_ft(0.5 * (horizon(l, t) + horizon(r, t)))
    };
    // Layer k (0 = top) is step s = nt − k, between horizon s (top) and s − 1 (base).
    let step_of = |k: usize| nt - k;

    let mut zcorn = Vec::with_capacity(8 * ni * nk);
    for k in 0..nk {
        let s = step_of(k);
        for t in [s, s - 1] {
            for _side_y in 0..2 {
                for i in 0..ni {
                    zcorn.push(edge(i, t));
                    zcorn.push(edge(i + 1, t));
                }
            }
        }
    }
    let ztop = zcorn.iter().cloned().fold(f64::INFINITY, f64::min);
    let zbot = zcorn.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut coord = Vec::with_capacity(6 * (ni + 1) * 2);
    let wy = opt.width_m * FT;
    for j in 0..2 {
        for e in 0..=ni {
            let x = (e as f64 * dx) * FT;
            let y = j as f64 * wy;
            coord.extend_from_slice(&[x, y, ztop - 10.0, x, y, zbot + 10.0]);
        }
    }
    let (mut actnum, mut poro, mut permx, mut permz, mut ntg) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut active_layers: Vec<Vec<usize>> = vec![Vec::new(); ni];
    let (mut pv, mut sum_phi, mut sum_ntg, mut n_act) = (0.0, 0.0, 0.0, 0usize);
    for k in 0..nk {
        let s = step_of(k);
        for (i, layers) in active_layers.iter_mut().enumerate() {
            let (h, phi, f) = rb.at(x0 + i, s);
            let active = h >= opt.min_thickness_m && phi > 0.0;
            let ks = kozeny(opt.sand_k_ref_md, phi, opt.phi_ref);
            let km = kozeny(opt.mud_k_ref_md, phi, opt.phi_ref);
            let k_mix = if active {
                (f * ks.ln() + (1.0 - f) * km.ln()).exp()
            } else {
                0.0
            };
            actnum.push(if active { 1.0 } else { 0.0 });
            poro.push(if active { phi } else { 0.0 });
            permx.push(k_mix);
            permz.push(k_mix * opt.kv_kh);
            ntg.push(if active { f.max(1e-3) } else { 0.0 });
            if active {
                layers.push(k);
                let vol_ft3 = h * FT * dx * FT * wy;
                pv += vol_ft3 * phi * f.max(1e-3) / 5.614_583;
                sum_phi += phi;
                sum_ntg += f;
                n_act += 1;
            }
        }
    }
    if n_act == 0 {
        return Err("window has no active cells".into());
    }
    let (inj, prod) = (0usize, ni - 1);
    if active_layers[inj].is_empty() || active_layers[prod].is_empty() {
        return Err("a well column has no active layers; choose another window".into());
    }
    let top_ft = ztop;
    let bot_ft = zbot;
    let summary = GridSummary {
        x_start_m: x0 as f64 * dx,
        x_end_m: (x0 + ni) as f64 * dx,
        nx: ni,
        nz: nk,
        active_cells: n_act,
        top_depth_ft: top_ft,
        bottom_depth_ft: bot_ft,
        pore_volume_rb: pv,
        mean_poro_active: sum_phi / n_act as f64,
        mean_ntg_active: sum_ntg / n_act as f64,
    };

    let mut d = String::new();
    let _ = writeln!(
        d,
        "-- Generated by strat-forward {} from a forward stratigraphic run.\n\
         -- Window {:.1}-{:.1} km, {} x 1 x {} cells (one layer per record step, youngest on top).\n\
         -- Rock model is a stand-in (Kozeny-Carman per class, log-linear mix); fluids are SPE1.\n",
        env!("CARGO_PKG_VERSION"),
        summary.x_start_m / 1000.0,
        summary.x_end_m / 1000.0,
        ni,
        nk
    );
    let _ = writeln!(
        d,
        "RUNSPEC\nTITLE\n  strat-forward section\nDIMENS\n  {ni} 1 {nk} /\nEQLDIMS\n/\nTABDIMS\n/\n\
         OIL\nGAS\nWATER\nDISGAS\nFIELD\nSTART\n  1 'JAN' 2030 /\nWELLDIMS\n  2 {nk} 1 2 /\n"
    );
    let _ = writeln!(d, "GRID\nSPECGRID\n  {ni} 1 {nk} 1 F /");
    write_array(&mut d, "COORD", &coord, 3);
    write_array(&mut d, "ZCORN", &zcorn, 3);
    write_array(&mut d, "ACTNUM", &actnum, 0);
    write_array(&mut d, "PORO", &poro, 4);
    write_array(&mut d, "PERMX", &permx, 4);
    write_array(&mut d, "PERMY", &permx, 4);
    write_array(&mut d, "PERMZ", &permz, 5);
    write_array(&mut d, "NTG", &ntg, 4);
    let _ = writeln!(d, "PINCH\n  {:.3} /\n", opt.min_thickness_m * FT);
    d.push_str(SPE1_PROPS);
    let owc = bot_ft + 100.0;
    let goc = top_ft - 100.0;
    let _ = writeln!(
        d,
        "SOLUTION\nEQUIL\n  {:.1} 4800 {:.1} 0 {:.1} 0 1 0 0 /\nRSVD\n  {:.1} 1.270\n  {:.1} 1.270 /\n",
        top_ft, owc, goc, goc, owc
    );
    let _ = writeln!(d, "SCHEDULE\nWELSPECS");
    let _ = writeln!(d, "  'PROD' 'G1' {} 1 {:.1} 'OIL' /", prod + 1, top_ft);
    let _ = writeln!(d, "  'INJ' 'G1' {} 1 {:.1} 'GAS' /\n/", inj + 1, top_ft);
    let _ = writeln!(d, "COMPDAT");
    for (well, col) in [("PROD", prod), ("INJ", inj)] {
        for k in &active_layers[col] {
            let _ = writeln!(
                d,
                "  '{well}' {} 1 {} {} 'OPEN' 1* 1* 0.5 /",
                col + 1,
                k + 1,
                k + 1
            );
        }
    }
    let _ = writeln!(d, "/");
    let _ = writeln!(
        d,
        "WCONPROD\n  'PROD' 'OPEN' 'ORAT' {:.1} 4* 1000 /\n/",
        opt.producer_oil_rate_stb_d
    );
    let _ = writeln!(
        d,
        "WCONINJE\n  'INJ' 'GAS' 'OPEN' 'RATE' {:.1} 1* 9014 /\n/",
        opt.injector_gas_rate_mscf_d
    );
    let _ = writeln!(d, "TSTEP\n  {}*30 /\nEND", opt.months);
    Ok(Deck { text: d, summary })
}

/// PROPS from the SPE1 benchmark (Odeh 1981), FIELD units.
const SPE1_PROPS: &str = "PROPS
PVTW
  4017.55 1.038 3.22E-6 0.318 0.0 /
ROCK
  14.7 3E-6 /
SWOF
0.12 0 1 0
0.18 4.64876033057851E-008 1 0
0.24 0.000000186 0.997 0
0.3 4.18388429752066E-007 0.98 0
0.36 7.43801652892562E-007 0.7 0
0.42 1.16219008264463E-006 0.35 0
0.48 1.67355371900826E-006 0.2 0
0.54 2.27789256198347E-006 0.09 0
0.6 2.97520661157025E-006 0.021 0
0.66 3.7654958677686E-006 0.01 0
0.72 4.64876033057851E-006 0.001 0
0.78 0.000005625 0.0001 0
0.84 6.69421487603306E-006 0 0
0.91 8.05914256198347E-006 0 0
1 0.00001 0 0 /
SGOF
0 0 1 0
0.001 0 1 0
0.02 0 0.997 0
0.05 0.005 0.980 0
0.12 0.025 0.700 0
0.2 0.075 0.350 0
0.25 0.125 0.200 0
0.3 0.190 0.090 0
0.4 0.410 0.021 0
0.45 0.60 0.010 0
0.5 0.72 0.001 0
0.6 0.87 0.0001 0
0.7 0.94 0.000 0
0.85 0.98 0.000 0
0.88 0.984 0.000 0 /
DENSITY
  53.66 64.49 0.0533 /
PVDG
14.700 166.666 0.008000
264.70 12.0930 0.009600
514.70 6.27400 0.011200
1014.7 3.19700 0.014000
2014.7 1.61400 0.018900
2514.7 1.29400 0.020800
3014.7 1.08000 0.022800
4014.7 0.81100 0.026800
5014.7 0.64900 0.030900
9014.7 0.38600 0.047000 /
PVTO
0.0010 14.7 1.0620 1.0400 /
0.0905 264.7 1.1500 0.9750 /
0.1800 514.7 1.2070 0.9100 /
0.3710 1014.7 1.2950 0.8300 /
0.6360 2014.7 1.4350 0.6950 /
0.7750 2514.7 1.5000 0.6410 /
0.9300 3014.7 1.5650 0.5940 /
1.2700 4014.7 1.6950 0.5100
  9014.7 1.5790 0.7400 /
1.6180 5014.7 1.8270 0.4490
  9014.7 1.7370 0.6310 /
/
";
