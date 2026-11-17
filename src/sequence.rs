//! Shoreline trajectory classes and the key sequence-stratigraphic turnarounds.
//!
//! The shoreline is the first sea-level crossing of the surface from the proximal end, interpolated between
//! cell centers. Each record step is classed by how the shoreline moved and what sea level did:
//! - forced regression: basinward while sea level falls
//! - normal regression: basinward while sea level rises or holds
//! - transgression: landward
//! - stationary: moved less than `min_move_m`
//!
//! Turnarounds are where the class changes. They are named for the surfaces they usually bound
//! (maximum flooding, maximum regression, onset and end of forced regression); whether a particular
//! turnaround deserves the name is an interpretation left to the reader.

use serde::Serialize;

use crate::model::{History, Model};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Trajectory {
    ForcedRegression,
    NormalRegression,
    Transgression,
    Stationary,
}

#[derive(Clone, Debug, Serialize)]
pub struct Step {
    pub step: usize,
    pub age_end_ma: f64,
    pub shoreline_x_m: Option<f64>,
    pub shoreline_dx_m: Option<f64>,
    pub sea_level_change_m: f64,
    pub class: Option<Trajectory>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Turnaround {
    /// The record step at whose end the change happens.
    pub step: usize,
    pub age_ma: f64,
    pub shoreline_x_m: Option<f64>,
    pub kind: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct Segment {
    pub class: Trajectory,
    pub first_step: usize,
    pub last_step: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct TrajectoryReport {
    pub min_move_m: f64,
    pub steps: Vec<Step>,
    pub segments: Vec<Segment>,
    pub turnarounds: Vec<Turnaround>,
}

/// Interpolated shoreline at record time `t`.
pub fn shoreline_at(model: &Model, h: &History, t: usize) -> Option<f64> {
    let nx = model.cfg.nx;
    let sl = h.sea_level_m[t];
    let eta = &h.eta_m[t * nx..(t + 1) * nx];
    let j = (0..nx).find(|&j| eta[j] < sl)?;
    let xc = |k: usize| (k as f64 + 0.5) * model.dx;
    if j == 0 {
        return Some(0.0);
    }
    let (a, b) = (eta[j - 1] - sl, eta[j] - sl); // a >= 0 > b
    Some(xc(j - 1) + model.dx * a / (a - b))
}

pub fn trajectory(model: &Model, h: &History) -> TrajectoryReport {
    let nt = model.cfg.n_records();
    let min_move_m = 0.05 * model.dx;
    let dt_ma = model.cfg.record_dt_yr / 1e6;
    let shore: Vec<Option<f64>> = (0..=nt).map(|t| shoreline_at(model, h, t)).collect();
    let mut steps = Vec::with_capacity(nt);
    for t in 1..=nt {
        let dsl = h.sea_level_m[t] - h.sea_level_m[t - 1];
        let dx = match (shore[t - 1], shore[t]) {
            (Some(a), Some(b)) => Some(b - a),
            _ => None,
        };
        let class = dx.map(|d| {
            if d.abs() < min_move_m {
                Trajectory::Stationary
            } else if d < 0.0 {
                Trajectory::Transgression
            } else if dsl < 0.0 {
                Trajectory::ForcedRegression
            } else {
                Trajectory::NormalRegression
            }
        });
        steps.push(Step {
            step: t,
            age_end_ma: model.cfg.start_age_ma - t as f64 * dt_ma,
            shoreline_x_m: shore[t],
            shoreline_dx_m: dx,
            sea_level_change_m: dsl,
            class,
        });
    }
    let mut segments: Vec<Segment> = Vec::new();
    for s in &steps {
        let Some(c) = s.class else { continue };
        match segments.last_mut() {
            Some(seg) if seg.class == c && seg.last_step + 1 == s.step => seg.last_step = s.step,
            _ => segments.push(Segment {
                class: c,
                first_step: s.step,
                last_step: s.step,
            }),
        }
    }
    // Turnarounds between moving segments; stationary segments are skipped over.
    let moving: Vec<&Segment> = segments
        .iter()
        .filter(|s| s.class != Trajectory::Stationary)
        .collect();
    let mut turnarounds = Vec::new();
    let regress = |c: Trajectory| {
        matches!(
            c,
            Trajectory::ForcedRegression | Trajectory::NormalRegression
        )
    };
    for pair in moving.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let at = a.last_step;
        // One turnaround can carry two names, e.g. transgression straight into forced regression is both
        // maximum flooding and the onset of forced regression.
        let mut kinds: Vec<&'static str> = Vec::new();
        if a.class == Trajectory::Transgression && regress(b.class) {
            kinds.push("maximum_flooding");
        }
        if regress(a.class) && b.class == Trajectory::Transgression {
            kinds.push("maximum_regression");
        }
        if a.class != Trajectory::ForcedRegression && b.class == Trajectory::ForcedRegression {
            kinds.push("onset_of_forced_regression");
        }
        if a.class == Trajectory::ForcedRegression && b.class != Trajectory::ForcedRegression {
            kinds.push("end_of_forced_regression");
        }
        for kind in kinds {
            turnarounds.push(Turnaround {
                step: at,
                age_ma: model.cfg.start_age_ma - at as f64 * dt_ma,
                shoreline_x_m: shore[at],
                kind,
            });
        }
    }
    TrajectoryReport {
        min_move_m,
        steps,
        segments,
        turnarounds,
    }
}
