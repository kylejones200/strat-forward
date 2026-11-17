//! The forward model: accommodation, conservative multi-class transport, compaction, and the event ledger.
//!
//! Three representations are kept apart:
//! - **State**: basement, surface, and the stack of units in each cell (`columns`). Compaction acts here.
//! - **Events**: the ordered ledger of what each cell did (`ledger`): deposit, erode, or exchange
//!   (one class in while another goes out), with the units each removal came from. Consecutive substeps
//!   doing the same thing to the same unit in the same record step are merged into one entry.
//! - **Stratigraphy**: the surviving units, each with a base contact (basement, conformable, diastem, or
//!   erosional with cut depth and the units it truncated). Derived from events, not from surface snapshots.
//!
//! A new unit starts at every tick of the internal layer clock (`layer_dt_yr`, which divides the record
//! interval), and also whenever the surface aggrades again after a cut deeper than `resolution_m`, so an
//! erosional contact survives even if deposition restores the elevation. Units are well mixed, so the layer
//! clock — never the record interval — sets the composition resolution transport sees.
//!
//! Transport is an explicit finite-volume diffusion of the surface. For class `i` at a face,
//! `q_i = −K_i(depth) · ∂η/∂x · w_i(source)`, where `w_i` is the class fraction in the active layer of the
//! upslope cell. Per-class flux leaving a cell in one substep is limited to what its active layer holds, so
//! solid volume is conserved exactly. Within a substep a cell's removals are applied before its deposits.
//! Between recompactions the surface moves by solid change / (1 − φ0); at each recompaction every column is
//! recompacted and the surface is reset to basement + compacted thickness.

use serde::Serialize;

use crate::compaction::bulk_thickness;
use crate::config::{Basement, Boundary, Config};

/// Solid thickness (m) below which a deposit is rounding noise and never starts a unit.
const DUST_M: f64 = 1e-12;

/// How a unit's base meets what is below it.
#[derive(Clone, Debug, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Contact {
    /// Directly on basement, with no measurable cut before it.
    Basement,
    /// On the unit below with no measurable cut and no missing record step.
    Conformable,
    /// No measurable cut, but whole record steps passed with no deposition.
    Diastem { missing_steps: usize },
    /// The surface was cut down by `cut_m` (decompacted) before this unit began.
    Erosional {
        cut_m: f64,
        truncated_units: Vec<u64>,
    },
}

/// One later removal from a unit, merged by record step and position.
#[derive(Clone, Debug, Serialize)]
pub struct Removal {
    pub step: usize,
    /// Solid removed per class.
    pub solid: Vec<f64>,
    /// True when a younger unit sat on top at the time: class exchange inside the active layer.
    pub beneath_younger: bool,
}

#[derive(Clone, Debug)]
pub struct Layer {
    /// Unique within the run.
    pub id: u64,
    /// Record step (1-based) during which this unit was laid down. 0 = initial sediment.
    pub step: usize,
    /// Tick of the layer clock during which this unit was laid down.
    pub slot: u64,
    /// Solid thickness per grain class still present (m).
    pub solid: Vec<f64>,
    /// Compacted thickness, refreshed at each recompaction.
    pub bulk: f64,
    /// Solid deposited per class, before any removal.
    pub deposited_class: Vec<f64>,
    pub t_first_yr: f64,
    pub t_last_yr: f64,
    pub base: Contact,
    /// Every later removal from this unit, in order.
    pub removals: Vec<Removal>,
    /// Origin tracer, per class: solid volume that had already been deposited in some unit and was eroded
    /// before arriving here (recycled), and Σ volume × time it entered the basin (yr from run start).
    /// Passive: transport never reads them. Initial sediment has entry time 0 and is not recycled.
    pub recycled: Vec<f64>,
    pub entry_moment: Vec<f64>,
    depth_moment: f64,
    deposited: f64,
}

impl Layer {
    pub fn total_solid(&self) -> f64 {
        self.solid.iter().sum()
    }
    /// Fraction of the surviving solid that is recycled sediment.
    pub fn recycled_fraction(&self) -> f64 {
        let tot = self.total_solid();
        if tot > 0.0 {
            self.recycled.iter().sum::<f64>() / tot
        } else {
            0.0
        }
    }
    /// Mean time (yr from run start) the surviving solid entered the basin.
    pub fn mean_entry_yr(&self) -> f64 {
        let tot = self.total_solid();
        if tot > 0.0 {
            self.entry_moment.iter().sum::<f64>() / tot
        } else {
            0.0
        }
    }
    /// Mean water depth at deposition, solid-weighted. Negative is above sea level.
    pub fn water_depth_m(&self) -> f64 {
        if self.deposited > 0.0 {
            self.depth_moment / self.deposited
        } else {
            0.0
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum OpKind {
    Deposit,
    Erode,
    /// Some classes arrive while others leave in the same substep.
    Exchange,
}

/// Solid taken from one unit by one ledger entry.
#[derive(Clone, Debug, Serialize)]
pub struct FromUnit {
    pub unit: u64,
    pub beneath_younger: bool,
    pub solid: Vec<f64>,
}

/// One ledger entry for one cell.
#[derive(Clone, Debug, Serialize)]
pub struct Op {
    pub step: usize,
    pub kind: OpKind,
    pub t_start_yr: f64,
    pub t_end_yr: f64,
    pub substeps: u32,
    /// Solid deposited per class, and the unit it went into.
    pub deposited: Vec<f64>,
    pub into_unit: Option<u64>,
    /// When this entry created `into_unit`: its base contact.
    pub created: Option<Contact>,
    pub removed: Vec<FromUnit>,
    /// Solid that crossed the cell without changing it: Σ min(in, out) · dt / dx.
    pub bypass_m: f64,
}

/// Mass budget per class, in solid m² (per unit width).
#[derive(Clone, Debug, serde::Serialize)]
pub struct Budget {
    pub class: String,
    pub initial_m2: f64,
    pub supplied_m2: f64,
    pub outflow_m2: f64,
    pub stored_m2: f64,
    /// initial + supplied − outflow − stored. Rounding only.
    pub error_m2: f64,
}

/// Per-cell bookkeeping for contacts: highest decompacted column since the current unit began, and the
/// lowest point since that high.
#[derive(Clone, Debug)]
struct Track {
    hi: f64,
    lo: f64,
    truncated: Vec<u64>,
    last_dep_step: usize,
}

pub struct Model {
    pub(crate) cfg: Config,
    pub(crate) dx: f64,
    pub(crate) basement: Vec<f64>,
    pub(crate) eta: Vec<f64>,
    pub(crate) columns: Vec<Vec<Layer>>,
    /// Ordered events per cell.
    pub(crate) ledger: Vec<Vec<Op>>,
    /// Initial units per cell (id, solid per class), the starting point for replaying the ledger.
    pub(crate) initial_units: Vec<Vec<(u64, Vec<f64>)>>,
    /// Solid per class per cell; kept equal to the layer sums.
    stock: Vec<Vec<f64>>,
    pub(crate) t_yr: f64,
    pub(crate) step: usize,
    /// Current tick of the layer clock.
    slot: u64,
    initial: Vec<f64>,
    supplied: Vec<f64>,
    outflow: Vec<f64>,
    substep_yr: f64,
    next_id: u64,
    track: Vec<Track>,
    // Scratch.
    flux: Vec<Vec<f64>>,
    weights: Vec<Vec<f64>>,
    /// Solid per class inside the active layer: the most a class can lose in one substep.
    active: Vec<Vec<f64>>,
    /// Decompacted thickness ever removed from each cell (substep resolution).
    eroded: Vec<f64>,
    /// Origin-tracer totals and per-substep scratch (see `trace`).
    pub(crate) tracer: TracerBudget,
    tr_eroded: Vec<Vec<[f64; 2]>>,
    /// Set once `run` has started; a model runs once.
    ran: bool,
    /// Flexural isostasy state, when configured.
    flex: Option<Flex>,
}

/// Plate state: column masses at the start (the reference load) and the current deflection.
struct Flex {
    rigidity: f64,
    rho_m: f64,
    rho_w: f64,
    rho_g: f64,
    response_time_yr: Option<f64>,
    /// Column mass per unit area at the start (kg/m²).
    mass0: Vec<f64>,
    /// Downward deflection since the start (m).
    w: Vec<f64>,
}

/// Origin-tracer totals per class (solid m² per unit width, entry moments in m²·yr). They close exactly:
/// stored + outflow entry moment = supplied entry moment (initial sediment enters at 0), and
/// stored + outflow recycled = recycled created by erosion.
#[derive(Clone, Debug, Default, Serialize)]
pub struct TracerBudget {
    pub supplied_entry_m2yr: Vec<f64>,
    pub outflow_entry_m2yr: Vec<f64>,
    pub outflow_recycled_m2: Vec<f64>,
    pub recycled_created_m2: Vec<f64>,
}

/// What the run recorded at every record time, oldest first (index 0 is the initial state).
pub struct History {
    pub time_yr: Vec<f64>,
    pub sea_level_m: Vec<f64>,
    /// Surface elevation, row-major `t * nx + x`.
    pub eta_m: Vec<f64>,
    /// Basement elevation, row-major.
    pub basement_m: Vec<f64>,
    /// Decompacted (deposited) sediment thickness above basement, row-major. Erosion lowers it; subsidence does not.
    pub deposited_m: Vec<f64>,
    /// Net decompacted change during step `t` (index t, t≥1), row-major. Row 0 is zero.
    pub net_m: Vec<f64>,
    /// Shoreline: first cell (from x=0) below sea level, or None if the section is all land.
    pub shoreline_x_m: Vec<Option<f64>>,
    /// Downward isostatic deflection of the basement since the start, row-major. Zero without isostasy.
    pub deflection_m: Vec<f64>,
}

impl Model {
    pub fn new(cfg: Config) -> crate::Result<Self> {
        cfg.validate()?;
        let nx = cfg.nx;
        let nc = cfg.classes.len();
        let dx = cfg.dx();
        let basement: Vec<f64> = match &cfg.initial.basement {
            Basement::Ramp { top_m, bottom_m } => (0..nx)
                .map(|j| top_m + (bottom_m - top_m) * (j as f64 + 0.5) / nx as f64)
                .collect(),
            Basement::Profile(p) => p.clone(),
        };
        let s0 = cfg.initial.sediment_solid_m.max(0.0);
        let mut columns = vec![Vec::new(); nx];
        let mut initial_units = vec![Vec::new(); nx];
        let mut stock = vec![vec![0.0; nc]; nx];
        let mut next_id = 1;
        if s0 > 0.0 {
            for j in 0..nx {
                let mut solid = vec![0.0; nc];
                solid[0] = s0;
                columns[j].push(Layer {
                    id: next_id,
                    step: 0,
                    slot: u64::MAX,
                    solid: solid.clone(),
                    bulk: 0.0,
                    deposited_class: solid.clone(),
                    t_first_yr: 0.0,
                    t_last_yr: 0.0,
                    base: Contact::Basement,
                    removals: Vec::new(),
                    recycled: vec![0.0; nc],
                    entry_moment: vec![0.0; nc],
                    depth_moment: 0.0,
                    deposited: s0,
                });
                initial_units[j].push((next_id, solid));
                next_id += 1;
                stock[j][0] = s0;
            }
        }
        let mut initial = vec![0.0; nc];
        initial[0] = s0 * nx as f64 * dx;

        // Explicit stability: the surface diffuses with at most K_i / (1 − φ0_i).
        let t = &cfg.transport;
        let d_max = cfg
            .classes
            .iter()
            .map(|c| {
                let k = (t.k_land_m2_per_yr * c.land_k_mult)
                    .max((t.k_marine_m2_per_yr + t.k_deep_m2_per_yr) * c.marine_k_mult);
                k / (1.0 - c.phi0)
            })
            .fold(0.0, f64::max);
        let substep_yr = if d_max > 0.0 {
            t.cfl * dx * dx / (2.0 * d_max)
        } else {
            cfg.record_dt_yr
        };

        let d0 = s0 / (1.0 - cfg.classes[0].phi0);
        let mut m = Model {
            dx,
            eta: basement.clone(),
            basement,
            columns,
            ledger: vec![Vec::new(); nx],
            initial_units,
            stock,
            t_yr: 0.0,
            step: 0,
            slot: 0,
            initial,
            supplied: vec![0.0; nc],
            outflow: vec![0.0; nc],
            substep_yr,
            next_id,
            track: vec![
                Track {
                    hi: d0,
                    lo: d0,
                    truncated: Vec::new(),
                    last_dep_step: 0
                };
                nx
            ],
            flux: vec![vec![0.0; nc]; nx + 1],
            weights: vec![vec![0.0; nc]; nx],
            active: vec![vec![0.0; nc]; nx],
            eroded: vec![0.0; nx],
            tracer: TracerBudget {
                supplied_entry_m2yr: vec![0.0; nc],
                outflow_entry_m2yr: vec![0.0; nc],
                outflow_recycled_m2: vec![0.0; nc],
                recycled_created_m2: vec![0.0; nc],
            },
            tr_eroded: vec![vec![[0.0; 2]; nc]; nx],
            ran: false,
            flex: None,
            cfg,
        };
        m.recompact();
        if let Some(iso) = m.cfg.isostasy.clone() {
            let mut f = Flex {
                rigidity: iso.rigidity(),
                rho_m: iso.mantle_density_kg_m3,
                rho_w: iso.water_density_kg_m3,
                rho_g: iso.grain_density_kg_m3,
                response_time_yr: iso.response_time_yr,
                mass0: Vec::new(),
                w: vec![0.0; nx],
            };
            let sl0 = m.cfg.sea_level.at(0.0);
            f.mass0 = (0..nx).map(|j| m.column_mass(&f, j, sl0)).collect();
            m.flex = Some(f);
        }
        Ok(m)
    }

    pub fn cfg(&self) -> &Config {
        &self.cfg
    }
    /// Cell width (m).
    pub fn dx(&self) -> f64 {
        self.dx
    }
    /// Basement elevation per cell (m).
    pub fn basement(&self) -> &[f64] {
        &self.basement
    }
    /// Surface elevation per cell (m).
    pub fn eta(&self) -> &[f64] {
        &self.eta
    }
    /// Surviving units per cell, oldest first.
    pub fn columns(&self) -> &[Vec<Layer>] {
        &self.columns
    }
    /// Ordered ledger entries per cell.
    pub fn ledger(&self) -> &[Vec<Op>] {
        &self.ledger
    }
    /// Initial units per cell (id, solid per class): the starting point for replaying the ledger.
    pub fn initial_units(&self) -> &[Vec<(u64, Vec<f64>)>] {
        &self.initial_units
    }
    /// Model time (yr from the start of the run).
    pub fn t_yr(&self) -> f64 {
        self.t_yr
    }
    /// Current record step (0 before the run).
    pub fn step(&self) -> usize {
        self.step
    }
    /// Origin-tracer totals.
    pub fn tracer(&self) -> &TracerBudget {
        &self.tracer
    }

    /// Decompacted thickness removed from each cell over the run, summed at substep resolution.
    pub fn eroded_m(&self) -> &[f64] {
        &self.eroded
    }

    pub fn substep_yr(&self) -> f64 {
        self.substep_yr
    }

    /// Downward isostatic deflection of the basement since the start (m), or None without isostasy.
    pub fn isostatic_deflection(&self) -> Option<&[f64]> {
        self.flex.as_ref().map(|f| f.w.as_slice())
    }

    /// Mass per unit area (kg/m²) above the basement in cell `j`: grains, pore water (pores assumed
    /// water-filled), and sea water above the bed when it is below sea level `sl`.
    fn column_mass(&self, f: &Flex, j: usize, sl: f64) -> f64 {
        let solid: f64 = self.stock[j].iter().sum();
        let bulk = (self.eta[j] - self.basement[j]).max(solid);
        f.rho_g * solid + f.rho_w * (bulk - solid) + f.rho_w * (sl - self.eta[j]).max(0.0)
    }

    /// Move the basement (and the surface with it) to the plate's deflection under the current load.
    ///
    /// Plate equilibrium is `D w'''' + ρ_m g w = g (m(w) − m0)`, with `m` the column mass above the
    /// deflected basement. Under water, deflecting further deepens the water by the same amount, so
    /// `m(w) = m(w_now) + ρ_w (w − w_now)` there; moving that term to the left gives the restoring force
    /// `(ρ_m − ρ_w) g` under water and `ρ_m g` on land. Exact while a cell stays on its side of sea level;
    /// a cell that crosses it is corrected on the next substep.
    fn flex_step(&mut self, dt: f64, sl: f64) {
        let Some(mut f) = self.flex.take() else {
            return;
        };
        let nx = self.cfg.nx;
        let mut restoring = vec![0.0; nx];
        let mut load = vec![0.0; nx];
        // `subsidence` is water-loaded tectonic subsidence (the backstripping convention): it already includes
        // the plate's response to water filling it, so that water is subtracted from the load whether or not
        // it is there. On land it is not, and the column rebounds by ρ_w T / ρ_m (the air-loaded equivalent).
        // Subtracting it on both sides keeps the equilibrium continuous at the shoreline; leaving it out on
        // land made it jump by ρ_w T / ρ_m there, and shoreline cells crossed back and forth every substep.
        let t = self.t_yr + dt;
        for j in 0..nx {
            let under = if self.eta[j] < sl { 1.0 } else { 0.0 };
            let tectonic = self.cfg.subsidence.rate((j as f64 + 0.5) / nx as f64) * t;
            restoring[j] = (f.rho_m - f.rho_w * under) * crate::isostasy::G;
            load[j] = crate::isostasy::G
                * (self.column_mass(&f, j, sl)
                    - f.mass0[j]
                    - f.rho_w * (under * f.w[j] + tectonic));
        }
        let target = crate::isostasy::deflection(self.dx, f.rigidity, &restoring, &load);
        let blend = f
            .response_time_yr
            .map_or(1.0, |tau| 1.0 - (-dt / tau).exp());
        for (j, tw) in target.iter().enumerate() {
            let dw = (tw - f.w[j]) * blend;
            f.w[j] += dw;
            self.basement[j] -= dw;
            self.eta[j] -= dw;
        }
        self.flex = Some(f);
    }

    fn diffusivity(&self, class: usize, water_depth: f64) -> f64 {
        let t = &self.cfg.transport;
        let c = &self.cfg.classes[class];
        let land = t.k_land_m2_per_yr * c.land_k_mult;
        let marine = |d: f64| {
            (t.k_marine_m2_per_yr * (-d.max(0.0) / t.depth_decay_m).exp() + t.k_deep_m2_per_yr)
                * c.marine_k_mult
        };
        if t.shoreline_width_m <= 0.0 {
            return if water_depth <= 0.0 {
                land
            } else {
                marine(water_depth)
            };
        }
        // Marine weight rises smoothly from 0 on land to 1 offshore, 0.5 at sea level.
        let s = 0.5 * (1.0 + (water_depth / t.shoreline_width_m).tanh());
        (1.0 - s) * land + s * marine(water_depth)
    }

    /// Per cell: `active` = solid of each class in the top `active_layer_m` (what may leave this substep), and
    /// `weights` = the class fractions transport sees.
    ///
    /// The fractions weight the column by depth below the surface with `exp(-z / active_layer_m)`. Under
    /// deposition this is exactly the composition of a well-mixed active layer of that thickness (each new
    /// deposit mixes in and the oldest material relaxes out exponentially). A sharp cut-off at
    /// `active_layer_m` instead makes the fractions a fixed window over the most recent deposits: once units
    /// are thinner than the window (a unit cap), they echo deposition one window-length earlier, and that
    /// delayed feedback oscillated, eroding up to 12 m and more on finer grids (examples/diag_cap_pattern.rs).
    fn active_weights(&mut self) {
        let active = self.cfg.transport.active_layer_m;
        for j in 0..self.cfg.nx {
            let a = &mut self.active[j];
            a.iter_mut().for_each(|v| *v = 0.0);
            let w = &mut self.weights[j];
            w.iter_mut().for_each(|v| *v = 0.0);
            let mut z = 0.0; // solid depth of this unit's top below the surface
            for layer in self.columns[j].iter().rev() {
                let tot = layer.total_solid();
                if tot <= 0.0 {
                    continue;
                }
                let in_active = (active - z).clamp(0.0, tot);
                let depth_weight = (-z / active).exp() - (-(z + tot) / active).exp();
                for ((ai, wi), si) in a.iter_mut().zip(w.iter_mut()).zip(&layer.solid) {
                    *ai += si * in_active / tot;
                    *wi += si / tot * depth_weight;
                }
                z += tot;
                // Below 30 e-folds the weight is under 1e-13 of the total.
                if z > 30.0 * active {
                    break;
                }
            }
            let sum: f64 = w.iter().sum();
            if sum > 0.0 {
                w.iter_mut().for_each(|v| *v /= sum);
            }
        }
    }

    fn substep(&mut self, dt: f64) {
        let nx = self.cfg.nx;
        let nc = self.cfg.classes.len();
        let dx = self.dx;
        let sl = self.cfg.sea_level.at(self.t_yr + 0.5 * dt);
        self.active_weights();

        // Face fluxes. flux[f] is the face on the proximal side of cell f; flux[nx] is the distal boundary.
        // Positive flux moves sediment distally (+x).
        for i in 0..nc {
            self.flux[0][i] = self.cfg.supply.flux_m2_per_yr * self.cfg.supply.fractions[i];
        }
        for f in 1..nx {
            let (a, b) = (f - 1, f);
            let slope = (self.eta[b] - self.eta[a]) / dx;
            let depth = sl - 0.5 * (self.eta[a] + self.eta[b]);
            let src = if slope < 0.0 { a } else { b };
            for i in 0..nc {
                self.flux[f][i] = -self.diffusivity(i, depth) * slope * self.weights[src][i];
            }
        }
        match self.cfg.transport.distal_boundary {
            Boundary::Closed => self.flux[nx].iter_mut().for_each(|q| *q = 0.0),
            Boundary::Open => {
                let slope = (self.eta[nx - 1] - self.eta[nx - 2]) / dx;
                let depth = sl - self.eta[nx - 1];
                for i in 0..nc {
                    self.flux[nx][i] = if slope < 0.0 {
                        -self.diffusivity(i, depth) * slope * self.weights[nx - 1][i]
                    } else {
                        0.0
                    };
                }
            }
        }

        // Limit what leaves each cell to what its active layer holds, so erosion never reaches below it.
        // Each face has one source, so scaling is conservative.
        for j in 0..nx {
            for i in 0..nc {
                let left = self.flux[j][i]; // out of j when negative
                let right = self.flux[j + 1][i]; // out of j when positive
                let out = ((-left).max(0.0) + right.max(0.0)) * dt / dx;
                let avail = self.active[j][i];
                if out > avail {
                    let scale = if out > 0.0 { avail / out } else { 0.0 };
                    if left < 0.0 && j > 0 {
                        self.flux[j][i] *= scale;
                    }
                    if right > 0.0 {
                        self.flux[j + 1][i] *= scale;
                    }
                }
            }
        }

        for i in 0..nc {
            self.supplied[i] += self.flux[0][i] * dt;
            self.outflow[i] += self.flux[nx][i] * dt;
        }

        let step = self.step;
        let tol = self.cfg.resolution_m;
        let t0 = self.t_yr;
        let mut ds = vec![0.0; nc];
        for row in self.tr_eroded.iter_mut() {
            row.iter_mut().for_each(|v| *v = [0.0; 2]);
        }
        for j in 0..nx {
            let depth = sl - self.eta[j];
            let mut bypass = 0.0;
            let mut any = false;
            for (i, d) in ds.iter_mut().enumerate() {
                let (l, r) = (self.flux[j][i], self.flux[j + 1][i]);
                let inflow = l.max(0.0) + (-r).max(0.0);
                let outflow = (-l).max(0.0) + r.max(0.0);
                bypass += inflow.min(outflow) * dt / dx;
                *d = (l - r) * dt / dx;
                any |= *d != 0.0;
            }
            if any || bypass > 0.0 {
                self.apply_cell(j, &ds, bypass, depth, step, tol, t0, dt);
            }
            let x_frac = (j as f64 + 0.5) / nx as f64;
            let sub = self.cfg.subsidence.rate(x_frac) * dt;
            self.basement[j] -= sub;
            self.eta[j] -= sub;
        }
        if self.flex.is_some() {
            self.flex_step(dt, self.cfg.sea_level.at(self.t_yr + dt));
        }
        self.trace(dt);
        self.t_yr += dt;
    }

    /// Carry the origin tracers through one substep, after the cells have been updated.
    ///
    /// Same convention as the ledger's bypass: in each cell and class, what leaves is first the material
    /// passing through, then what was eroded from the bed; what is deposited is the part of the inflow that
    /// stayed. Outflow is well mixed and splits between the two faces in proportion to the flux. With X_j
    /// the tracer content leaving cell j per unit time, b_j = min(in, out)/in the share of inflow that passes
    /// through, and ρ, λ the shares of a cell's outflow going distal and proximal:
    ///   X_j = b_j (ρ_{j−1} X_{j−1} + λ_{j+1} X_{j+1} + S_j) + C_j
    /// with S the supply and C the content of eroded bed material. Each face carries flow one way, so the
    /// flow graph has no cycles and this tridiagonal system is nonsingular. Content is conserved exactly
    /// and every value is a mix of its sources.
    fn trace(&mut self, dt: f64) {
        let nx = self.cfg.nx;
        let nc = self.cfg.classes.len();
        let dx = self.dx;
        let t_mid = self.t_yr + 0.5 * dt;
        let (mut a, mut c, mut d) = (vec![0.0; nx], vec![0.0; nx], vec![[0.0f64; 2]; nx]);
        let (mut rho, mut lam, mut byp, mut dep) =
            (vec![0.0; nx], vec![0.0; nx], vec![0.0; nx], vec![0.0; nx]);
        for i in 0..nc {
            for j in 0..nx {
                let (l, r) = (self.flux[j][i], self.flux[j + 1][i]);
                let inflow = l.max(0.0) + (-r).max(0.0);
                let out = (-l).max(0.0) + r.max(0.0);
                rho[j] = if out > 0.0 { r.max(0.0) / out } else { 0.0 };
                lam[j] = if out > 0.0 { (-l).max(0.0) / out } else { 0.0 };
                byp[j] = if inflow > 0.0 {
                    inflow.min(out) / inflow
                } else {
                    0.0
                };
                dep[j] = if inflow > out { 1.0 - byp[j] } else { 0.0 };
            }
            // Supply enters through face 0: not recycled, entering now.
            let q0 = self.flux[0][i].max(0.0);
            self.tracer.supplied_entry_m2yr[i] += q0 * t_mid * dt;
            let supply = [0.0, q0 * t_mid];
            for j in 0..nx {
                a[j] = if j > 0 { -byp[j] * rho[j - 1] } else { 0.0 };
                c[j] = if j + 1 < nx {
                    -byp[j] * lam[j + 1]
                } else {
                    0.0
                };
                // Eroded bed material, as a rate: all recycled; entry moment from the units it left.
                let e = self.tr_eroded[j][i];
                let mut rhs = [e[0] * dx / dt, e[1] * dx / dt];
                if j == 0 {
                    rhs[0] += byp[0] * supply[0];
                    rhs[1] += byp[0] * supply[1];
                }
                d[j] = rhs;
            }
            let x = thomas(&a, &c, &d);
            for j in 0..nx {
                if dep[j] <= 0.0 {
                    continue;
                }
                let mut inc = [0.0; 2];
                for k in 0..2 {
                    if j > 0 {
                        inc[k] += rho[j - 1] * x[j - 1][k];
                    }
                    if j + 1 < nx {
                        inc[k] += lam[j + 1] * x[j + 1][k];
                    }
                    if j == 0 {
                        inc[k] += supply[k];
                    }
                }
                let top = self.columns[j].last_mut().expect("a deposit made a unit");
                top.recycled[i] += dep[j] * inc[0] * dt / dx;
                top.entry_moment[i] += dep[j] * inc[1] * dt / dx;
            }
            let last = nx - 1;
            self.tracer.outflow_recycled_m2[i] += rho[last] * x[last][0] * dt;
            self.tracer.outflow_entry_m2yr[i] += rho[last] * x[last][1] * dt;
        }
    }

    /// Apply one substep's per-class changes to cell `j`: removals first, then deposits; update contacts and ledger.
    #[allow(clippy::too_many_arguments)]
    fn apply_cell(
        &mut self,
        j: usize,
        ds: &[f64],
        bypass: f64,
        depth: f64,
        step: usize,
        tol: f64,
        t0: f64,
        dt: f64,
    ) {
        let nc = ds.len();
        let phi0: Vec<f64> = self.cfg.classes.iter().map(|c| c.phi0).collect();
        let mut d_eta = 0.0;
        let mut from: Vec<FromUnit> = Vec::new();

        // Removals.
        for i in 0..nc {
            if ds[i] < 0.0 {
                let (removed, rec_removed, entry_removed) = erode(
                    &mut self.columns[j],
                    i,
                    -ds[i],
                    step,
                    &mut from,
                    &mut self.track[j].truncated,
                );
                // Eroded material is recycled from now on, whatever it was before.
                self.tr_eroded[j][i][0] += removed;
                self.tr_eroded[j][i][1] += entry_removed;
                self.tracer.recycled_created_m2[i] += (removed - rec_removed) * self.dx;
                self.stock[j][i] = (self.stock[j][i] - removed).max(0.0);
                d_eta -= removed / (1.0 - phi0[i]);
                self.eroded[j] += removed / (1.0 - phi0[i]);
            }
        }
        let d_mid = decompacted_with(&self.stock[j], &phi0);
        {
            let tr = &mut self.track[j];
            tr.lo = tr.lo.min(d_mid);
        }

        // Deposits.
        let dep_dec: f64 = (0..nc)
            .filter(|&i| ds[i] > 0.0)
            .map(|i| ds[i] / (1.0 - phi0[i]))
            .sum();
        let net = dep_dec + d_eta;
        let mut created: Option<Contact> = None;
        // A rounding-level deposit (supply and outflow nearly canceling) joins the top unit instead of
        // starting one: such a unit would hold ~1e-18 m, be dropped later with no ledger entry, and take with
        // it any erosional contact owed to the next real unit.
        let dep_solid: f64 = ds.iter().filter(|v| **v > 0.0).sum();
        if dep_dec > 0.0 && dep_solid < DUST_M && !self.columns[j].is_empty() {
            let top = self.columns[j].last_mut().expect("unit");
            d_eta += deposit_into(top, &mut self.stock[j], ds, &phi0, depth);
        } else if dep_dec > 0.0 {
            let tr = &self.track[j];
            let cut = tr.hi - tr.lo;
            let col = &self.columns[j];
            let slot = self.slot;
            let cap = self.cfg.max_unit_solid_m.unwrap_or(f64::INFINITY);
            let need_new = col.is_empty()
                || col
                    .last()
                    .is_none_or(|l| l.slot != slot || l.total_solid() >= cap)
                || (net > 0.0 && cut > tol);
            if need_new {
                let contact = if cut > tol {
                    Contact::Erosional {
                        cut_m: cut,
                        truncated_units: tr.truncated.clone(),
                    }
                } else if col.is_empty() {
                    Contact::Basement
                } else if step > tr.last_dep_step + 1 {
                    Contact::Diastem {
                        missing_steps: step - tr.last_dep_step - 1,
                    }
                } else {
                    Contact::Conformable
                };
                let id = self.next_id;
                self.next_id += 1;
                self.columns[j].push(Layer {
                    id,
                    step,
                    slot,
                    solid: vec![0.0; nc],
                    bulk: 0.0,
                    deposited_class: vec![0.0; nc],
                    t_first_yr: t0,
                    t_last_yr: t0 + dt,
                    base: contact.clone(),
                    removals: Vec::new(),
                    recycled: vec![0.0; nc],
                    entry_moment: vec![0.0; nc],
                    depth_moment: 0.0,
                    deposited: 0.0,
                });
                created = Some(contact);
                self.track[j].truncated.clear();
            }
            let top = self.columns[j].last_mut().expect("unit");
            d_eta += deposit_into(top, &mut self.stock[j], ds, &phi0, depth);
            top.t_last_yr = t0 + dt;
            self.track[j].last_dep_step = step;
        }
        self.eta[j] += d_eta;

        let d_after = decompacted_with(&self.stock[j], &phi0);
        let tr = &mut self.track[j];
        if created.is_some() || d_after >= tr.hi {
            tr.hi = d_after;
            tr.lo = d_after;
        } else {
            tr.lo = tr.lo.min(d_after);
        }

        // Ledger.
        let has_dep = dep_dec > 0.0;
        let has_rem = !from.is_empty();
        let kind = match (has_dep, has_rem) {
            (true, true) => OpKind::Exchange,
            (false, true) => OpKind::Erode,
            _ => OpKind::Deposit, // deposit, or pure bypass with no change
        };
        let into = if has_dep {
            self.columns[j].last().map(|l| l.id)
        } else {
            None
        };
        let deposited: Vec<f64> = ds.iter().map(|v| v.max(0.0)).collect();
        // A rounding-level change (supply and outflow canceling to ~1e-17 m in a bypass cell) joins the
        // previous entry of the same step whatever its kind. As separate entries these flipped between
        // deposit and erode every substep and made up 24% of the reference ledger. Replay is unaffected: it adds
        // each entry's removals and deposits and does not read the kind.
        let removed_solid: f64 = from.iter().flat_map(|f| f.solid.iter()).sum();
        let dust = dep_solid < DUST_M && removed_solid < DUST_M;
        let ops = &mut self.ledger[j];
        let merge = created.is_none()
            && ops.last().is_some_and(|o| {
                o.step == step
                    && if dust {
                        into.is_none() || o.into_unit.is_none() || o.into_unit == into
                    } else {
                        o.kind == kind && o.into_unit == into
                    }
            });
        if merge {
            let o = ops.last_mut().expect("op");
            if o.into_unit.is_none() {
                o.into_unit = into;
            }
            o.t_end_yr = t0 + dt;
            o.substeps += 1;
            o.bypass_m += bypass;
            for (sum, d) in o.deposited.iter_mut().zip(&deposited) {
                *sum += d;
            }
            for f in from {
                merge_from(&mut o.removed, f);
            }
        } else {
            ops.push(Op {
                step,
                kind,
                t_start_yr: t0,
                t_end_yr: t0 + dt,
                substeps: 1,
                deposited,
                into_unit: into,
                created,
                removed: from,
                bypass_m: bypass,
            });
        }
    }

    /// Recompute every layer's compacted thickness and reset the surface.
    pub fn recompact(&mut self) {
        let classes = &self.cfg.classes;
        for (j, col) in self.columns.iter_mut().enumerate() {
            let mut z = 0.0;
            for layer in col.iter_mut().rev() {
                let tot = layer.total_solid();
                if tot <= 0.0 {
                    layer.bulk = 0.0;
                    continue;
                }
                let (mut phi0, mut c) = (0.0, 0.0);
                for (s, cls) in layer.solid.iter().zip(classes) {
                    phi0 += s / tot * cls.phi0;
                    c += s / tot * cls.compaction_length_m;
                }
                layer.bulk = bulk_thickness(tot, z, phi0, c);
                z += layer.bulk;
            }
            self.eta[j] = self.basement[j] + z;
        }
    }

    /// Decompacted (as-deposited) thickness of a solid vector.
    pub fn decompacted(&self, solid: &[f64]) -> f64 {
        solid
            .iter()
            .zip(&self.cfg.classes)
            .map(|(s, c)| s / (1.0 - c.phi0))
            .sum()
    }

    fn shoreline(&self, sl: f64) -> Option<f64> {
        (0..self.cfg.nx)
            .find(|&j| self.eta[j] < sl)
            .map(|j| (j as f64 + 0.5) * self.dx)
    }

    /// Run to the end and return the record. A model runs once; a second call is an error.
    pub fn run(&mut self) -> crate::Result<History> {
        if self.ran {
            return Err(crate::Error::State(
                "Model::run called twice; build a new Model to run again".into(),
            ));
        }
        self.ran = true;
        let nx = self.cfg.nx;
        let nt = self.cfg.n_records();
        let layer_dt = self.cfg.layer_dt();
        let mut next_tick = layer_dt;
        let mut h = History {
            time_yr: Vec::with_capacity(nt + 1),
            sea_level_m: Vec::with_capacity(nt + 1),
            eta_m: Vec::with_capacity((nt + 1) * nx),
            basement_m: Vec::with_capacity((nt + 1) * nx),
            deposited_m: Vec::with_capacity((nt + 1) * nx),
            net_m: vec![0.0; nx],
            shoreline_x_m: Vec::with_capacity(nt + 1),
            deflection_m: Vec::with_capacity((nt + 1) * nx),
        };
        self.snapshot(&mut h);
        for rec in 1..=nt {
            self.step = rec;
            let before: Vec<f64> = (0..nx).map(|j| self.decompacted(&self.stock[j])).collect();
            let t_end = rec as f64 * self.cfg.record_dt_yr;
            while self.t_yr < t_end - 1e-9 {
                let dt = self
                    .substep_yr
                    .min(t_end - self.t_yr)
                    .min(next_tick - self.t_yr);
                self.substep(dt);
                if self.t_yr >= next_tick - 1e-9 {
                    // Tick: recompact and start new units from here on. The record end is always a tick.
                    self.t_yr = next_tick;
                    self.recompact();
                    self.slot += 1;
                    next_tick += layer_dt;
                }
            }
            self.t_yr = t_end;
            for (stock, b) in self.stock.iter().zip(&before) {
                h.net_m.push(self.decompacted(stock) - b);
            }
            self.snapshot(&mut h);
        }
        Ok(h)
    }

    fn snapshot(&self, h: &mut History) {
        let sl = self.cfg.sea_level.at(self.t_yr);
        h.time_yr.push(self.t_yr);
        h.sea_level_m.push(sl);
        h.eta_m.extend_from_slice(&self.eta);
        h.basement_m.extend_from_slice(&self.basement);
        for j in 0..self.cfg.nx {
            h.deposited_m.push(self.decompacted(&self.stock[j]));
        }
        h.shoreline_x_m.push(self.shoreline(sl));
        match &self.flex {
            Some(f) => h.deflection_m.extend_from_slice(&f.w),
            None => h.deflection_m.extend(std::iter::repeat_n(0.0, self.cfg.nx)),
        }
    }

    pub fn budget(&self) -> Vec<Budget> {
        let nc = self.cfg.classes.len();
        (0..nc)
            .map(|i| {
                let stored: f64 = self
                    .columns
                    .iter()
                    .map(|col| col.iter().map(|l| l.solid[i]).sum::<f64>())
                    .sum::<f64>()
                    * self.dx;
                Budget {
                    class: self.cfg.classes[i].name.clone(),
                    initial_m2: self.initial[i],
                    supplied_m2: self.supplied[i],
                    outflow_m2: self.outflow[i],
                    stored_m2: stored,
                    error_m2: self.initial[i] + self.supplied[i] - self.outflow[i] - stored,
                }
            })
            .collect()
    }
}

fn decompacted_with(solid: &[f64], phi0: &[f64]) -> f64 {
    solid.iter().zip(phi0).map(|(s, p)| s / (1.0 - p)).sum()
}

fn merge_from(list: &mut Vec<FromUnit>, f: FromUnit) {
    if let Some(e) = list
        .iter_mut()
        .find(|e| e.unit == f.unit && e.beneath_younger == f.beneath_younger)
    {
        for (a, b) in e.solid.iter_mut().zip(&f.solid) {
            *a += b;
        }
    } else {
        list.push(f);
    }
}

/// Add this substep's positive per-class changes `ds` to `top` and the cell's `stock`. Returns the
/// decompacted thickness added.
fn deposit_into(top: &mut Layer, stock: &mut [f64], ds: &[f64], phi0: &[f64], depth: f64) -> f64 {
    let mut d_eta = 0.0;
    for (i, &d) in ds.iter().enumerate() {
        if d > 0.0 {
            top.solid[i] += d;
            top.deposited_class[i] += d;
            top.deposited += d;
            top.depth_moment += d * depth;
            stock[i] += d;
            d_eta += d / (1.0 - phi0[i]);
        }
    }
    d_eta
}

/// Remove up to `amount` of one class from the top down, recording which units lost it. Returns the solid
/// removed and the recycled volume and entry moment it carried out of the units.
fn erode(
    col: &mut Vec<Layer>,
    class: usize,
    amount: f64,
    step: usize,
    from: &mut Vec<FromUnit>,
    truncated: &mut Vec<u64>,
) -> (f64, f64, f64) {
    let (mut rec_out, mut entry_out) = (0.0, 0.0);
    let mut need = amount;
    let top = col.len().saturating_sub(1);
    for (k, layer) in col.iter_mut().enumerate().rev() {
        if need <= 0.0 {
            break;
        }
        let take = layer.solid[class].min(need);
        if take <= 0.0 {
            continue;
        }
        let share = take / layer.solid[class];
        let (dr, de) = (
            layer.recycled[class] * share,
            layer.entry_moment[class] * share,
        );
        layer.recycled[class] -= dr;
        layer.entry_moment[class] -= de;
        rec_out += dr;
        entry_out += de;
        layer.solid[class] -= take;
        need -= take;
        let beneath = k != top;
        let mut solid = vec![0.0; layer.solid.len()];
        solid[class] = take;
        match layer.removals.last_mut() {
            Some(r) if r.step == step && r.beneath_younger == beneath => r.solid[class] += take,
            _ => layer.removals.push(Removal {
                step,
                solid: solid.clone(),
                beneath_younger: beneath,
            }),
        }
        merge_from(
            from,
            FromUnit {
                unit: layer.id,
                beneath_younger: beneath,
                solid,
            },
        );
        if !truncated.contains(&layer.id) {
            truncated.push(layer.id);
        }
    }
    col.retain(|l| l.total_solid() > 1e-15);
    (amount - need, rec_out, entry_out)
}

/// Solve a tridiagonal system with unit diagonal for two right-hand sides at once. `a` is the sub- and `c`
/// the super-diagonal. No pivoting: the tracer matrix is column diagonally dominant.
fn thomas(a: &[f64], c: &[f64], d: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let n = d.len();
    let mut cp = vec![0.0; n];
    let mut dp = vec![[0.0; 2]; n];
    let mut denom = 1.0;
    cp[0] = c[0] / denom;
    dp[0] = [d[0][0] / denom, d[0][1] / denom];
    for k in 1..n {
        denom = 1.0 - a[k] * cp[k - 1];
        cp[k] = c[k] / denom;
        for r in 0..2 {
            dp[k][r] = (d[k][r] - a[k] * dp[k - 1][r]) / denom;
        }
    }
    let mut x = vec![[0.0; 2]; n];
    x[n - 1] = dp[n - 1];
    for k in (0..n - 1).rev() {
        for r in 0..2 {
            x[k][r] = dp[k][r] - cp[k] * x[k + 1][r];
        }
    }
    x
}
