//! Portable outputs: JSON datasets and one SVG figure (final section over a Wheeler diagram).

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::error::Error;
use crate::model::{Budget, History, Model};
use crate::record::{self, Well, Wheeler};

/// Elevation is decompacted sediment thickness above basement: the frame without later subsidence or
/// compaction, so the suffix minimum over time is exactly the preserved top of each surface.
#[derive(Serialize)]
pub struct SurfaceGridOut {
    pub length_m: f64,
    pub resolution_m: f64,
    pub nx: usize,
    pub elevation_m: Vec<f64>,
    pub time_ma: Vec<f64>,
    pub age_ma: Vec<f64>,
}

pub fn surface_grid(model: &Model, h: &History) -> SurfaceGridOut {
    let cfg = &model.cfg;
    SurfaceGridOut {
        length_m: cfg.length_m,
        resolution_m: cfg.resolution_m,
        nx: cfg.nx,
        elevation_m: h.deposited_m.clone(),
        time_ma: h.time_yr.iter().map(|t| t / 1e6).collect(),
        age_ma: h
            .time_yr
            .iter()
            .map(|t| cfg.start_age_ma - t / 1e6)
            .collect(),
    }
}

#[derive(Serialize)]
struct Summary<'a> {
    engine: &'static str,
    version: &'static str,
    config: &'a crate::Config,
    substep_yr: f64,
    budget: Vec<Budget>,
    max_relative_mass_error: f64,
    surface_reconstruction_gap_m: f64,
    gap_breakdown: GapBreakdown,
    ledger_replay_max_diff_m: f64,
    ledger_entries: usize,
    surviving_units: usize,
    contacts: record::ContactCounts,
    time_yr: &'a [f64],
    sea_level_m: &'a [f64],
    shoreline_x_m: &'a [Option<f64>],
    wheeler_counts: &'a record::Counts,
    /// Origin-tracer totals (see `Model::tracer`) and the recycled share of all surviving solid.
    tracer: &'a crate::model::TracerBudget,
    recycled_fraction: f64,
    /// Present only with isostasy, so outputs without it are unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    isostasy: Option<IsostasySummary>,
}

/// The plate and how far it moved.
#[derive(Serialize)]
struct IsostasySummary {
    rigidity_nm: f64,
    /// Flexural parameter under water, (4 D / ((rho_m - rho_w) g))^(1/4): how far a load's effect reaches.
    flexural_parameter_m: f64,
    max_deflection_m: f64,
    /// Per record time, row-major `t * nx + x`, like the other grids.
    deflection_m: Vec<f64>,
}

/// Share of all surviving solid that is recycled sediment (eroded from a unit and redeposited).
pub fn recycled_fraction(model: &Model) -> f64 {
    let (mut rec, mut tot) = (0.0, 0.0);
    for l in model.columns.iter().flatten() {
        rec += l.recycled.iter().sum::<f64>();
        tot += l.total_solid();
    }
    if tot > 0.0 {
        rec / tot
    } else {
        0.0
    }
}

/// Largest per-step difference (m, decompacted) between what a surface-only reconstruction (the suffix
/// minimum over the sampled surface grid) and the event record say a
/// step preserved. Non-zero when cut-and-fill happens between samples or when grain classes exchange
/// within the active layer.
pub fn surface_reconstruction_gap(model: &Model, h: &History) -> f64 {
    let (nx, nt) = (model.cfg.nx, model.cfg.n_records());
    let w = record::wheeler(model, h);
    let mut worst = 0.0_f64;
    for x in 0..nx {
        let surfaces: Vec<f64> = (0..=nt).map(|t| h.deposited_m[t * nx + x]).collect();
        let pres = record::suffix_min(&surfaces);
        for t in 1..=nt {
            let from_surfaces = (pres[t] - pres[t - 1]).max(0.0);
            worst = worst.max((from_surfaces - w.preserved_m[(t - 1) * nx + x]).abs());
        }
    }
    worst
}

/// Where the surface/record gap comes from, over every (step, cell) whose gap exceeds `resolution_m`.
#[derive(Clone, Debug, Serialize)]
pub struct GapBreakdown {
    pub threshold_m: f64,
    pub cell_steps_over_threshold: usize,
    pub cell_steps_total: usize,
    pub dominated_by_exchange_beneath_younger: usize,
    pub dominated_by_later_top_erosion: usize,
    pub dominated_by_cut_into_older_within_step: usize,
    pub max_identity_residual_m: f64,
    pub worst: Option<record::GapTerms>,
}

pub fn gap_breakdown(model: &Model, h: &History) -> GapBreakdown {
    let g = record::gap_terms(model, h);
    let tol = model.cfg.resolution_m;
    let gap = |t: &record::GapTerms| (t.surface_m.max(0.0) - t.record_m).abs();
    let big: Vec<&record::GapTerms> = g.iter().filter(|t| gap(t) > tol).collect();
    let count = |f: &dyn Fn(&record::GapTerms) -> bool| big.iter().filter(|t| f(t)).count();
    GapBreakdown {
        threshold_m: tol,
        cell_steps_over_threshold: big.len(),
        cell_steps_total: g.len(),
        dominated_by_exchange_beneath_younger: count(&|t| {
            t.removed_later_beneath_m >= t.removed_later_top_m.max(t.cut_into_older_m)
        }),
        dominated_by_later_top_erosion: count(&|t| {
            t.removed_later_top_m > t.removed_later_beneath_m.max(t.cut_into_older_m)
        }),
        dominated_by_cut_into_older_within_step: count(&|t| {
            t.cut_into_older_m > t.removed_later_beneath_m.max(t.removed_later_top_m)
        }),
        max_identity_residual_m: g
            .iter()
            .map(|t| t.identity_residual_m.abs())
            .fold(0.0, f64::max),
        worst: g.iter().max_by(|a, b| gap(a).total_cmp(&gap(b))).cloned(),
    }
}

#[derive(Serialize)]
struct WheelerOut<'a> {
    #[serde(flatten)]
    wheeler: &'a Wheeler,
    /// Gross decompacted deposition and erosion per (step, cell) from the ledger. Both can be non-zero in a
    /// step whose net change is small: cut-and-fill the state codes alone do not show.
    gross_deposited_m: Vec<f64>,
    gross_eroded_m: Vec<f64>,
}

pub fn max_relative_mass_error(budget: &[Budget]) -> f64 {
    budget
        .iter()
        .map(|b| b.error_m2.abs() / (b.initial_m2 + b.supplied_m2).abs().max(1e-30))
        .fold(0.0, f64::max)
}

pub fn default_wells(model: &Model) -> Vec<Well> {
    let nx = model.cfg.nx;
    [0.15, 0.35, 0.55, 0.75]
        .iter()
        .enumerate()
        .map(|(k, f)| {
            record::well(
                model,
                &format!("W-{}", k + 1),
                ((nx as f64 * f) as usize).min(nx - 1),
            )
        })
        .collect()
}

/// Write every standard output into `dir` (created if missing). JSON is streamed to disk, never built as one
/// string in memory. Returns the paths written, in order.
pub fn write_all(dir: &Path, model: &Model, h: &History) -> crate::Result<Vec<PathBuf>> {
    std::fs::create_dir_all(dir).map_err(|e| Error::io(dir, e))?;
    let wheeler = record::wheeler(model, h);
    let wells = default_wells(model);
    let budget = model.budget();
    let strat = record::stratigraphy(model);
    let replay = record::replay(model)?;
    let (gross_dep, gross_ero) = record::gross(model);
    let summary = Summary {
        engine: env!("CARGO_PKG_NAME"),
        version: env!("CARGO_PKG_VERSION"),
        config: &model.cfg,
        substep_yr: model.substep_yr(),
        max_relative_mass_error: max_relative_mass_error(&budget),
        surface_reconstruction_gap_m: surface_reconstruction_gap(model, h),
        gap_breakdown: gap_breakdown(model, h),
        ledger_replay_max_diff_m: replay,
        ledger_entries: model.ledger.iter().map(|o| o.len()).sum(),
        surviving_units: model.columns.iter().map(|c| c.len()).sum(),
        contacts: strat.contacts.clone(),
        budget,
        time_yr: &h.time_yr,
        sea_level_m: &h.sea_level_m,
        shoreline_x_m: &h.shoreline_x_m,
        wheeler_counts: &wheeler.counts,
        tracer: &model.tracer,
        recycled_fraction: recycled_fraction(model),
        isostasy: model.cfg.isostasy.as_ref().map(|iso| {
            let d = iso.rigidity();
            let restoring =
                (iso.mantle_density_kg_m3 - iso.water_density_kg_m3) * crate::isostasy::G;
            IsostasySummary {
                rigidity_nm: d,
                flexural_parameter_m: crate::isostasy::flexural_parameter(d, restoring),
                max_deflection_m: h.deflection_m.iter().cloned().fold(0.0, f64::max),
                deflection_m: h.deflection_m.clone(),
            }
        }),
    };
    let mut written = Vec::new();
    written.push(write_json(dir, "summary.json", &summary, true)?);
    written.push(write_json(
        dir,
        "surface_grid.json",
        &surface_grid(model, h),
        true,
    )?);
    let wheeler_out = WheelerOut {
        wheeler: &wheeler,
        gross_deposited_m: gross_dep,
        gross_eroded_m: gross_ero,
    };
    written.push(write_json(dir, "wheeler.json", &wheeler_out, true)?);
    written.push(write_json(dir, "stratigraphy.json", &strat, true)?);
    written.push(write_json(
        dir,
        "trajectory.json",
        &crate::sequence::trajectory(model, h),
        true,
    )?);
    // The ledger is large; written compact.
    written.push(write_json(dir, "ledger.json", &model.ledger, false)?);
    written.push(write_json(dir, "wells.json", &wells, true)?);
    written.push(write_json(
        dir,
        "hiatus_surfaces.json",
        &crate::record::hiatus_surfaces(model, h),
        true,
    )?);
    written.push(write_with(dir, "section.svg", |w| {
        w.write_all(svg(model, h, &wheeler, &wells).as_bytes())
    })?);
    Ok(written)
}

/// Stream `value` as JSON to `dir/name`.
pub fn write_json<T: serde::Serialize + ?Sized>(
    dir: &Path,
    name: &str,
    value: &T,
    pretty: bool,
) -> crate::Result<PathBuf> {
    write_with(dir, name, |w| {
        let r = if pretty {
            serde_json::to_writer_pretty(&mut *w, value)
        } else {
            serde_json::to_writer(&mut *w, value)
        };
        r.map_err(std::io::Error::other)
    })
}

/// Create `dir/name` and fill it through a buffered writer; flushes and reports any I/O error with the path.
pub fn write_with(
    dir: &Path,
    name: &str,
    fill: impl FnOnce(&mut std::io::BufWriter<std::fs::File>) -> std::io::Result<()>,
) -> crate::Result<PathBuf> {
    let p = dir.join(name);
    let f = std::fs::File::create(&p).map_err(|e| Error::io(&p, e))?;
    let mut w = std::io::BufWriter::new(f);
    fill(&mut w)
        .and_then(|_| w.flush())
        .map_err(|e| Error::io(&p, e))?;
    Ok(p)
}

fn mix(a: (f64, f64, f64), b: (f64, f64, f64), t: f64) -> String {
    let t = t.clamp(0.0, 1.0);
    let c = |x: f64, y: f64| (x + (y - x) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", c(a.0, b.0), c(a.1, b.1), c(a.2, b.2))
}

const MUD: (f64, f64, f64) = (96.0, 104.0, 112.0);
const SAND: (f64, f64, f64) = (232.0, 196.0, 84.0);

fn facies_color(f: &str) -> &'static str {
    match f {
        "coastal plain" => "#7fae6a",
        "shoreface" => "#e8c454",
        "shelf" => "#6fa3c8",
        _ => "#45607a",
    }
}

/// Final cross section (colored by sand fraction, timelines every 5 steps) over a Wheeler diagram.
pub fn svg(model: &Model, h: &History, w: &Wheeler, wells: &[Well]) -> String {
    let cfg = &model.cfg;
    let nx = cfg.nx;
    let nt = cfg.n_records();
    let (width, left, right) = (1200.0, 70.0, 30.0);
    let pw = width - left - right;
    let (sec_top, sec_h) = (40.0, 380.0);
    let (wh_top, wh_h) = (sec_top + sec_h + 70.0, 300.0);
    let height = wh_top + wh_h + 60.0;
    let dx = model.dx;
    let xs = |x: f64| left + x / cfg.length_m * pw;

    let last = nt;
    let sl_now = h.sea_level_m[last];
    let zmax = model.eta.iter().cloned().fold(sl_now, f64::max) + 10.0;
    let zmin = model.basement.iter().cloned().fold(f64::INFINITY, f64::min) - 10.0;
    let ys = |z: f64| sec_top + (zmax - z) / (zmax - zmin) * sec_h;

    let mut s = String::new();
    let _ = write!(
        s,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="Helvetica, Arial, sans-serif" font-size="11">
<rect width="100%" height="100%" fill="#fbfaf7"/>
<text x="{left}" y="22" font-size="14" font-weight="bold">Final dip section — color = sand fraction, black = timelines every 5 steps, red = erosional contacts ({:.2} Ma)</text>
<rect x="{left}" y="{sec_top}" width="{pw}" height="{sec_h}" fill="#dfe9f1"/>
"##,
        5.0 * cfg.record_dt_yr / 1e6
    );
    // Water above the final surface.
    let _ = writeln!(
        s,
        r##"<rect x="{left}" y="{}" width="{pw}" height="{}" fill="#bcd4e6"/>"##,
        ys(sl_now),
        (sec_top + sec_h - ys(sl_now)).max(0.0)
    );
    // Layers.
    let mut timeline: Vec<Vec<(f64, f64)>> = vec![Vec::new(); nt / 5 + 1];
    let mut contact_segs: Vec<String> = Vec::new();
    for (j, col) in model.columns.iter().enumerate() {
        let x0 = xs(j as f64 * dx);
        let x1 = xs((j + 1) as f64 * dx);
        let mut z = model.basement[j];
        let mut k = 0;
        for (tl, line) in timeline.iter_mut().enumerate() {
            let step_lim = tl * 5;
            // Top of all layers with step <= step_lim.
            let mut zz = model.basement[j];
            for l in col {
                if l.step <= step_lim {
                    zz += l.bulk;
                }
            }
            line.push(((x0 + x1) / 2.0, ys(zz)));
        }
        for l in col {
            let top = z + l.bulk;
            if l.bulk > 1e-3 {
                let _ = writeln!(
                    s,
                    r##"<rect x="{x0:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{}"/>"##,
                    ys(top),
                    x1 - x0 + 0.3,
                    (ys(z) - ys(top)).max(0.1),
                    mix(MUD, SAND, record::sand_fraction(l))
                );
            }
            if matches!(l.base, crate::model::Contact::Erosional { .. }) && l.bulk > 1e-3 {
                contact_segs.push(format!(
                    r##"<line x1="{x0:.2}" x2="{x1:.2}" y1="{y:.2}" y2="{y:.2}" stroke="#d7261e" stroke-width="1.6"/>"##,
                    y = ys(z)
                ));
            }
            z = top;
            k += 1;
        }
        let _ = k;
    }
    for line in timeline.iter().skip(1) {
        let pts: Vec<String> = line.iter().map(|(x, y)| format!("{x:.1},{y:.1}")).collect();
        let _ = writeln!(
            s,
            r##"<polyline points="{}" fill="none" stroke="#222" stroke-width="0.6" opacity="0.7"/>"##,
            pts.join(" ")
        );
    }
    for seg in &contact_segs {
        s.push_str(seg);
        s.push('\n');
    }
    let base: Vec<String> = (0..nx)
        .map(|j| {
            format!(
                "{:.1},{:.1}",
                xs((j as f64 + 0.5) * dx),
                ys(model.basement[j])
            )
        })
        .collect();
    let _ = writeln!(
        s,
        r##"<polyline points="{}" fill="none" stroke="#5a3b2e" stroke-width="1.6"/>"##,
        base.join(" ")
    );
    let _ = writeln!(
        s,
        r##"<line x1="{left}" x2="{}" y1="{y:.1}" y2="{y:.1}" stroke="#1f6fb2" stroke-dasharray="5 3"/><text x="{}" y="{:.1}" fill="#1f6fb2" text-anchor="end">sea level {sl_now:.1} m</text>"##,
        left + pw,
        left + pw - 4.0,
        ys(sl_now) - 4.0,
        y = ys(sl_now)
    );
    for wl in wells {
        let x = xs(wl.x_m);
        let _ = writeln!(
            s,
            r##"<line x1="{x:.1}" x2="{x:.1}" y1="{:.1}" y2="{:.1}" stroke="#c0392b" stroke-width="1.2"/><text x="{x:.1}" y="{:.1}" fill="#c0392b" text-anchor="middle">{}</text>"##,
            ys(wl.surface_elevation_m) - 8.0,
            ys(wl.basement_elevation_m),
            ys(wl.surface_elevation_m) - 11.0,
            wl.name
        );
    }
    // Depth axis.
    let step_z = nice_step((zmax - zmin) / 6.0);
    let mut z = (zmin / step_z).ceil() * step_z;
    while z <= zmax {
        let _ = writeln!(
            s,
            r##"<line x1="{}" x2="{left}" y1="{y:.1}" y2="{y:.1}" stroke="#333"/><text x="{}" y="{:.1}" text-anchor="end">{z:.0}</text>"##,
            left - 4.0,
            left - 6.0,
            ys(z) + 4.0,
            y = ys(z)
        );
        z += step_z;
    }
    let _ = writeln!(
        s,
        r##"<text x="14" y="{:.0}" transform="rotate(-90 14 {:.0})" text-anchor="middle">elevation (m)</text>"##,
        sec_top + sec_h / 2.0,
        sec_top + sec_h / 2.0
    );

    // Wheeler.
    let age0 = cfg.start_age_ma;
    let age1 = cfg.start_age_ma - nt as f64 * cfg.record_dt_yr / 1e6;
    let ya = |a: f64| wh_top + (age0 - a) / (age0 - age1) * wh_h;
    let _ = writeln!(
        s,
        r##"<text x="{left}" y="{:.0}" font-size="14" font-weight="bold">Wheeler diagram — preserved deposits by facies; gray hiatus, red erosion, pink vacuity (deposited then removed)</text>"##,
        wh_top - 14.0
    );
    for t in 0..nt {
        let y0 = ya(w.age_start_ma[t]);
        let y1 = ya(w.age_end_ma[t]);
        for x in 0..nx {
            let at = t * nx + x;
            let fill = match w.state[at] {
                1 => facies_color(w.facies[at].unwrap_or("shelf")),
                2 => "#e4e4e4",
                3 => "#c8553d",
                _ => "#f2b8c6",
            };
            let _ = writeln!(
                s,
                r##"<rect x="{:.2}" y="{y0:.2}" width="{:.2}" height="{:.2}" fill="{fill}"/>"##,
                xs(x as f64 * dx),
                xs(dx) - left + 0.3,
                y1 - y0 + 0.3
            );
        }
    }
    // Shoreline trajectory on the Wheeler.
    let pts: Vec<String> = h
        .shoreline_x_m
        .iter()
        .zip(&h.time_yr)
        .filter_map(|(sx, t)| sx.map(|x| format!("{:.1},{:.1}", xs(x), ya(age0 - t / 1e6))))
        .collect();
    let _ = writeln!(
        s,
        r##"<polyline points="{}" fill="none" stroke="#111" stroke-width="1.4"/>"##,
        pts.join(" ")
    );
    let step_a = nice_step((age0 - age1) / 6.0);
    let mut a = (age1 / step_a).ceil() * step_a;
    while a <= age0 + 1e-9 {
        let _ = writeln!(
            s,
            r##"<line x1="{}" x2="{left}" y1="{y:.1}" y2="{y:.1}" stroke="#333"/><text x="{}" y="{:.1}" text-anchor="end">{a:.1}</text>"##,
            left - 4.0,
            left - 6.0,
            ya(a) + 4.0,
            y = ya(a)
        );
        a += step_a;
    }
    let _ = writeln!(
        s,
        r##"<text x="14" y="{:.0}" transform="rotate(-90 14 {:.0})" text-anchor="middle">age (Ma)</text>"##,
        wh_top + wh_h / 2.0,
        wh_top + wh_h / 2.0
    );
    let step_x = nice_step(cfg.length_m / 1000.0 / 8.0);
    let mut xk = 0.0;
    while xk <= cfg.length_m / 1000.0 + 1e-9 {
        let x = xs(xk * 1000.0);
        let _ = writeln!(
            s,
            r##"<line x1="{x:.1}" x2="{x:.1}" y1="{}" y2="{}" stroke="#333"/><text x="{x:.1}" y="{:.0}" text-anchor="middle">{xk:.0}</text>"##,
            wh_top + wh_h,
            wh_top + wh_h + 4.0,
            wh_top + wh_h + 16.0
        );
        xk += step_x;
    }
    let _ = writeln!(
        s,
        r##"<text x="{:.0}" y="{:.0}" text-anchor="middle">distance basinward (km) — black line: shoreline trajectory</text>"##,
        left + pw / 2.0,
        wh_top + wh_h + 34.0
    );
    let mut lx = left;
    for (name, c) in [
        ("coastal plain", "#7fae6a"),
        ("shoreface", "#e8c454"),
        ("shelf", "#6fa3c8"),
        ("slope-basin", "#45607a"),
        ("hiatus", "#e4e4e4"),
        ("erosion", "#c8553d"),
        ("vacuity", "#f2b8c6"),
    ] {
        let _ = writeln!(
            s,
            r##"<rect x="{lx}" y="{:.0}" width="12" height="12" fill="{c}" stroke="#999"/><text x="{}" y="{:.0}">{name}</text>"##,
            height - 20.0,
            lx + 16.0,
            height - 10.0
        );
        lx += 110.0;
    }
    s.push_str("</svg>\n");
    s
}

fn nice_step(raw: f64) -> f64 {
    let p = 10f64.powf(raw.abs().max(1e-12).log10().floor());
    let m = raw / p;
    p * if m < 1.5 {
        1.0
    } else if m < 3.5 {
        2.0
    } else if m < 7.5 {
        5.0
    } else {
        10.0
    }
}
