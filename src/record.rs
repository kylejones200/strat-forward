//! Readouts of a finished run: chronostratigraphic occupancy (Wheeler), synthetic wells, facies.
//!
//! Occupancy uses four states (deposition, hiatus, erosion, vacuity) and is read from the event
//! record rather than inferred from surfaces.

use serde::Serialize;

use crate::model::{History, Layer, Model};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ChronoState {
    Deposition,
    Hiatus,
    Erosion,
    Vacuity,
}

impl ChronoState {
    pub fn code(self) -> u8 {
        match self {
            ChronoState::Deposition => 1,
            ChronoState::Hiatus => 2,
            ChronoState::Erosion => 3,
            ChronoState::Vacuity => 4,
        }
    }
    pub fn decision(self) -> &'static str {
        match self {
            ChronoState::Deposition => "rise_preserved",
            ChronoState::Hiatus => "within_resolution",
            ChronoState::Erosion => "surface_lowered",
            ChronoState::Vacuity => "rise_removed",
        }
    }
}

pub fn facies(water_depth_m: f64) -> &'static str {
    if water_depth_m <= 0.0 {
        "coastal plain"
    } else if water_depth_m < 15.0 {
        "shoreface"
    } else if water_depth_m < 80.0 {
        "shelf"
    } else {
        "slope-basin"
    }
}

pub fn lithology(sand_fraction: f64) -> &'static str {
    if sand_fraction >= 0.6 {
        "sandstone"
    } else if sand_fraction >= 0.3 {
        "heterolithic"
    } else {
        "mudstone"
    }
}

/// Fraction of class 0 (by convention the coarsest class) in a layer's solids.
pub fn sand_fraction(layer: &Layer) -> f64 {
    let tot = layer.total_solid();
    if tot > 0.0 {
        layer.solid[0] / tot
    } else {
        0.0
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Wheeler {
    pub nx: usize,
    /// Number of record steps. Row `t` (0-based) is step `t + 1`.
    pub nt: usize,
    pub x_m: Vec<f64>,
    /// Age at the start and end of each step (Ma).
    pub age_start_ma: Vec<f64>,
    pub age_end_ma: Vec<f64>,
    /// State code per (step, x), row-major: 1 deposition, 2 hiatus, 3 erosion, 4 vacuity.
    pub state: Vec<u8>,
    /// Preserved decompacted thickness of each step's deposit at the end of the run.
    pub preserved_m: Vec<f64>,
    /// Net decompacted change during the step.
    pub net_m: Vec<f64>,
    /// Facies of the preserved deposit, or null.
    pub facies: Vec<Option<&'static str>>,
    pub counts: Counts,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Counts {
    pub deposition: usize,
    pub hiatus: usize,
    pub erosion: usize,
    pub vacuity: usize,
}

pub fn classify(state_net: f64, preserved: f64, tol: f64) -> ChronoState {
    if state_net > tol {
        if preserved > tol {
            ChronoState::Deposition
        } else {
            ChronoState::Vacuity
        }
    } else if state_net < -tol {
        ChronoState::Erosion
    } else {
        ChronoState::Hiatus
    }
}

pub fn wheeler(model: &Model, h: &History) -> Wheeler {
    let cfg = &model.cfg;
    let nx = cfg.nx;
    let nt = cfg.n_records();
    let tol = cfg.resolution_m;
    let mut preserved = vec![0.0; nt * nx];
    let mut fac = vec![None; nt * nx];
    for (x, col) in model.columns.iter().enumerate() {
        for layer in col {
            if layer.step >= 1 {
                let at = (layer.step - 1) * nx + x;
                preserved[at] += model.decompacted(&layer.solid);
                fac[at] = Some(facies(layer.water_depth_m()));
            }
        }
    }
    let mut state = Vec::with_capacity(nt * nx);
    let mut net = Vec::with_capacity(nt * nx);
    let mut counts = Counts::default();
    for t in 0..nt {
        for x in 0..nx {
            let n = h.net_m[(t + 1) * nx + x];
            let s = classify(n, preserved[t * nx + x], tol);
            match s {
                ChronoState::Deposition => counts.deposition += 1,
                ChronoState::Hiatus => counts.hiatus += 1,
                ChronoState::Erosion => counts.erosion += 1,
                ChronoState::Vacuity => counts.vacuity += 1,
            }
            if s != ChronoState::Deposition {
                fac[t * nx + x] = None;
            }
            state.push(s.code());
            net.push(n);
        }
    }
    let dt_ma = cfg.record_dt_yr / 1e6;
    Wheeler {
        nx,
        nt,
        x_m: (0..nx).map(|j| (j as f64 + 0.5) * model.dx).collect(),
        age_start_ma: (0..nt)
            .map(|t| cfg.start_age_ma - t as f64 * dt_ma)
            .collect(),
        age_end_ma: (0..nt)
            .map(|t| cfg.start_age_ma - (t + 1) as f64 * dt_ma)
            .collect(),
        state,
        preserved_m: preserved,
        net_m: net,
        facies: fac,
        counts,
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct WellInterval {
    pub step: usize,
    pub age_top_ma: f64,
    pub age_base_ma: f64,
    pub depth_top_m: f64,
    pub depth_base_m: f64,
    pub thickness_m: f64,
    pub porosity: f64,
    pub sand_fraction: f64,
    pub water_depth_at_deposition_m: f64,
    pub facies: &'static str,
    pub lithology: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct Well {
    pub name: String,
    pub x_m: f64,
    pub cell: usize,
    pub surface_elevation_m: f64,
    pub basement_elevation_m: f64,
    /// Youngest first.
    pub intervals: Vec<WellInterval>,
}

/// A synthetic well at cell `x`, intervals youngest first, depths below the final surface.
pub fn well(model: &Model, name: &str, x: usize) -> Well {
    let cfg = &model.cfg;
    let dt_ma = cfg.record_dt_yr / 1e6;
    let mut z = 0.0;
    let mut intervals = Vec::new();
    for layer in model.columns[x].iter().rev() {
        let s = layer.total_solid();
        if s <= 0.0 {
            continue;
        }
        let sf = sand_fraction(layer);
        let wd = layer.water_depth_m();
        let step = layer.step;
        intervals.push(WellInterval {
            step,
            age_top_ma: cfg.start_age_ma - step as f64 * dt_ma,
            age_base_ma: cfg.start_age_ma - (step.max(1) - 1) as f64 * dt_ma,
            depth_top_m: z,
            depth_base_m: z + layer.bulk,
            thickness_m: layer.bulk,
            porosity: crate::compaction::mean_porosity(s, layer.bulk),
            sand_fraction: sf,
            water_depth_at_deposition_m: wd,
            facies: facies(wd),
            lithology: lithology(sf),
        });
        z += layer.bulk;
    }
    Well {
        name: name.to_string(),
        x_m: (x as f64 + 0.5) * model.dx,
        cell: x,
        surface_elevation_m: model.eta[x],
        basement_elevation_m: model.basement[x],
        intervals,
    }
}

// ---------------------------------------------------------------------------------------------------------
// Readouts built from the event ledger.

use std::collections::HashMap;

use crate::model::Contact;

/// Replay every cell's ledger from its initial units and compare with the final column.
/// Returns the largest per-class solid difference (m). Errors if the ledger names a unit it never created,
/// or the surviving unit order differs.
pub fn replay(model: &Model) -> crate::Result<f64> {
    replay_inner(model).map_err(crate::Error::Invariant)
}

fn replay_inner(model: &Model) -> Result<f64, String> {
    let mut worst = 0.0_f64;
    for (x, ops) in model.ledger.iter().enumerate() {
        let mut order: Vec<u64> = Vec::new();
        let mut solids: HashMap<u64, Vec<f64>> = HashMap::new();
        for (id, s) in &model.initial_units[x] {
            order.push(*id);
            solids.insert(*id, s.clone());
        }
        for (k, op) in ops.iter().enumerate() {
            if op.created.is_some() {
                let id = op
                    .into_unit
                    .ok_or_else(|| format!("cell {x} op {k}: created without a unit"))?;
                order.push(id);
                solids.insert(id, vec![0.0; op.deposited.len()]);
            }
            for f in &op.removed {
                let s = solids.get_mut(&f.unit).ok_or_else(|| {
                    format!("cell {x} op {k}: removal from unknown unit {}", f.unit)
                })?;
                for (a, b) in s.iter_mut().zip(&f.solid) {
                    *a -= b;
                }
            }
            if let Some(id) = op.into_unit {
                let s = solids
                    .get_mut(&id)
                    .ok_or_else(|| format!("cell {x} op {k}: deposit into unknown unit {id}"))?;
                for (a, b) in s.iter_mut().zip(&op.deposited) {
                    *a += b;
                }
            }
        }
        let survivors: Vec<u64> = order
            .into_iter()
            .filter(|id| solids[id].iter().sum::<f64>() > 1e-12)
            .collect();
        let actual: Vec<u64> = model.columns[x]
            .iter()
            .filter(|l| l.total_solid() > 1e-12)
            .map(|l| l.id)
            .collect();
        if survivors != actual {
            return Err(format!(
                "cell {x}: replayed units {:?} != model units {:?}",
                survivors, actual
            ));
        }
        for l in &model.columns[x] {
            if let Some(s) = solids.get(&l.id) {
                for (a, b) in s.iter().zip(&l.solid) {
                    worst = worst.max((a - b).abs());
                }
            }
        }
    }
    Ok(worst)
}

/// For one (step, cell): how a surface-only reconstruction and the event record differ, term by term.
/// All values are decompacted thickness (m). With
///   surface = pres[t] − pres[t−1]   (suffix minimum of sampled surfaces, before clamping at zero)
///   record  = deposited − removed_same_step − removed_later_top − removed_later_beneath
/// the identity  surface − record = removed_later − later_lowering_t + later_lowering_t_minus_1 − cut_into_older
/// holds exactly; `identity_residual` is what is left of it numerically.
#[derive(Clone, Debug, Serialize)]
pub struct GapTerms {
    pub x: usize,
    pub step: usize,
    pub surface_m: f64,
    pub record_m: f64,
    /// Deposited into this step's units during the step.
    pub deposited_m: f64,
    /// Removed from this step's units during the same step.
    pub removed_same_step_m: f64,
    /// Removed from older units during this step: cut-and-fill hidden inside the step.
    pub cut_into_older_m: f64,
    /// Removed from this step's units in later steps while they were on top.
    pub removed_later_top_m: f64,
    /// Removed from this step's units in later steps while a younger unit lay above (class exchange).
    pub removed_later_beneath_m: f64,
    /// How far later sampled surfaces fall below this step's sampled surface.
    pub later_lowering_m: f64,
    /// Same for the previous step's surface.
    pub later_lowering_prev_m: f64,
    pub identity_residual_m: f64,
}

/// Gap terms for every (step, cell), row-major `(step − 1) * nx + x`.
pub fn gap_terms(model: &Model, h: &History) -> Vec<GapTerms> {
    let nx = model.cfg.nx;
    let nt = model.cfg.n_records();
    let mut out = Vec::with_capacity(nt * nx);
    let mut per_cell: Vec<Vec<GapTerms>> = Vec::with_capacity(nx);
    for x in 0..nx {
        let mut unit_step: HashMap<u64, usize> = model.initial_units[x]
            .iter()
            .map(|(id, _)| (*id, 0))
            .collect();
        for op in &model.ledger[x] {
            if let (Some(_), Some(id)) = (&op.created, op.into_unit) {
                unit_step.insert(id, op.step);
            }
        }
        let (mut g, mut e_in, mut e_old, mut l_top, mut l_bel) = (
            vec![0.0; nt + 1],
            vec![0.0; nt + 1],
            vec![0.0; nt + 1],
            vec![0.0; nt + 1],
            vec![0.0; nt + 1],
        );
        for op in &model.ledger[x] {
            g[op.step] += model.decompacted(&op.deposited);
            for f in &op.removed {
                let us = unit_step[&f.unit];
                let d = model.decompacted(&f.solid);
                if us == op.step {
                    e_in[op.step] += d;
                } else if us < op.step {
                    e_old[op.step] += d;
                    if us >= 1 {
                        if f.beneath_younger {
                            l_bel[us] += d;
                        } else {
                            l_top[us] += d;
                        }
                    }
                }
            }
        }
        let d: Vec<f64> = (0..=nt).map(|t| h.deposited_m[t * nx + x]).collect();
        let pres = suffix_min(&d);
        let mut col = Vec::with_capacity(nt);
        for t in 1..=nt {
            let surface = pres[t] - pres[t - 1];
            let record = g[t] - e_in[t] - l_top[t] - l_bel[t];
            let v = d[t] - pres[t];
            let w = d[t - 1] - pres[t - 1];
            let rhs = l_top[t] + l_bel[t] - v + w - e_old[t];
            col.push(GapTerms {
                x,
                step: t,
                surface_m: surface,
                record_m: record,
                deposited_m: g[t],
                removed_same_step_m: e_in[t],
                cut_into_older_m: e_old[t],
                removed_later_top_m: l_top[t],
                removed_later_beneath_m: l_bel[t],
                later_lowering_m: v,
                later_lowering_prev_m: w,
                identity_residual_m: (surface - record) - rhs,
            });
        }
        per_cell.push(col);
    }
    for t in 0..nt {
        for cell in &per_cell {
            out.push(cell[t].clone());
        }
    }
    out
}

#[derive(Clone, Debug, Serialize)]
pub struct UnitOut {
    pub id: u64,
    pub step: usize,
    pub age_start_ma: f64,
    pub age_end_ma: f64,
    pub t_first_yr: f64,
    pub t_last_yr: f64,
    pub base: Contact,
    pub deposited_m: f64,
    pub preserved_m: f64,
    pub preserved_fraction: f64,
    pub compacted_m: f64,
    pub sand_fraction: f64,
    pub water_depth_at_deposition_m: f64,
    pub facies: &'static str,
    /// Later events that removed part of this unit: (step, decompacted m, beneath a younger unit).
    pub modified_by: Vec<(usize, f64, bool)>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ContactCounts {
    pub basement: usize,
    pub conformable: usize,
    pub diastem: usize,
    pub erosional: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Stratigraphy {
    pub nx: usize,
    pub nt: usize,
    /// Surviving units per cell, oldest first.
    pub columns: Vec<Vec<UnitOut>>,
    pub contacts: ContactCounts,
    /// Final elevation of the top of all surviving units laid down by step t, row-major `t * nx + x`.
    /// A preserved horizon, not a snapshot: it is where the record of time t lies now.
    pub preserved_horizon_m: Vec<f64>,
}

pub fn stratigraphy(model: &Model) -> Stratigraphy {
    let cfg = &model.cfg;
    let nx = cfg.nx;
    let nt = cfg.n_records();
    let dt_ma = cfg.record_dt_yr / 1e6;
    let mut contacts = ContactCounts::default();
    let mut columns = Vec::with_capacity(nx);
    let mut horizon = vec![0.0; (nt + 1) * nx];
    for (x, col) in model.columns.iter().enumerate() {
        let mut units = Vec::with_capacity(col.len());
        for l in col {
            match l.base {
                Contact::Basement => contacts.basement += 1,
                Contact::Conformable => contacts.conformable += 1,
                Contact::Diastem { .. } => contacts.diastem += 1,
                Contact::Erosional { .. } => contacts.erosional += 1,
            }
            let dep = model.decompacted(&l.deposited_class);
            let pre = model.decompacted(&l.solid);
            let sf = sand_fraction(l);
            let wd = l.water_depth_m();
            units.push(UnitOut {
                id: l.id,
                step: l.step,
                age_start_ma: cfg.start_age_ma - (l.step.max(1) - 1) as f64 * dt_ma,
                age_end_ma: cfg.start_age_ma - l.step as f64 * dt_ma,
                t_first_yr: l.t_first_yr,
                t_last_yr: l.t_last_yr,
                base: l.base.clone(),
                deposited_m: dep,
                preserved_m: pre,
                preserved_fraction: if dep > 0.0 { pre / dep } else { 0.0 },
                compacted_m: l.bulk,
                sand_fraction: sf,
                water_depth_at_deposition_m: wd,
                facies: facies(wd),
                modified_by: l
                    .removals
                    .iter()
                    .map(|r| (r.step, model.decompacted(&r.solid), r.beneath_younger))
                    .collect(),
            });
        }
        for t in 0..=nt {
            let top: f64 = col.iter().filter(|l| l.step <= t).map(|l| l.bulk).sum();
            horizon[t * nx + x] = model.basement[x] + top;
        }
        columns.push(units);
    }
    Stratigraphy {
        nx,
        nt,
        columns,
        contacts,
        preserved_horizon_m: horizon,
    }
}

/// Gross deposition and erosion per (step, cell), decompacted m, row-major `(step − 1) * nx + x`.
pub fn gross(model: &Model) -> (Vec<f64>, Vec<f64>) {
    let nx = model.cfg.nx;
    let nt = model.cfg.n_records();
    let (mut dep, mut ero) = (vec![0.0; nt * nx], vec![0.0; nt * nx]);
    for (x, ops) in model.ledger.iter().enumerate() {
        for op in ops {
            if op.step == 0 {
                continue;
            }
            let at = (op.step - 1) * nx + x;
            dep[at] += model.decompacted(&op.deposited);
            ero[at] += op
                .removed
                .iter()
                .map(|f| model.decompacted(&f.solid))
                .sum::<f64>();
        }
    }
    (dep, ero)
}

/// A run of hiatus record steps at one cell, with preserved deposition directly below and above it.
#[derive(Clone, Debug, Serialize)]
pub struct HiatusSurface {
    pub x: usize,
    /// Last deposition step below the gap.
    pub below_step: usize,
    /// First deposition step above the gap.
    pub above_step: usize,
    pub missing_steps: usize,
}

/// Diastems read from the Wheeler states, so sub-resolution trickle deposition does not hide them.
pub fn hiatus_surfaces(model: &Model, h: &History) -> Vec<HiatusSurface> {
    let w = wheeler(model, h);
    let (nx, nt) = (w.nx, w.nt);
    let mut out = Vec::new();
    for x in 0..nx {
        let mut last_dep: Option<usize> = None;
        let mut gap = 0usize;
        for t in 0..nt {
            match w.state[t * nx + x] {
                1 => {
                    if let (Some(b), true) = (last_dep, gap > 0) {
                        out.push(HiatusSurface {
                            x,
                            below_step: b + 1,
                            above_step: t + 1,
                            missing_steps: gap,
                        });
                    }
                    last_dep = Some(t);
                    gap = 0;
                }
                2 => gap += 1,
                _ => {
                    // Erosion or vacuity: the gap is an erosional surface, recorded by the units' contacts.
                    last_dep = None;
                    gap = 0;
                }
            }
        }
    }
    out
}

/// Preserved top of each sampled surface: the minimum of it and every later surface. This is what a
/// surface-only reconstruction sees as the record.
pub(crate) fn suffix_min(surfaces: &[f64]) -> Vec<f64> {
    let mut pres = surfaces.to_vec();
    for t in (0..pres.len().saturating_sub(1)).rev() {
        pres[t] = pres[t].min(pres[t + 1]);
    }
    pres
}
