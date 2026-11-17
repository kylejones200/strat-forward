//! Run configuration. Units: meters, years. Fluxes are solid volume per unit width (m²/yr).

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub length_m: f64,
    pub nx: usize,
    pub duration_yr: f64,
    /// Interval between recorded surfaces. Each record step is one layer in the stratigraphic record.
    pub record_dt_yr: f64,
    /// Age of the start of the run, in Ma. Ages decrease as the model runs.
    pub start_age_ma: f64,
    /// Smallest decompacted thickness change (m) treated as deposition or erosion. Also the smallest cut
    /// that starts a new unit with an erosional base.
    pub resolution_m: f64,
    /// The model's internal clock: a new unit starts, and every column is recompacted, at each multiple.
    /// Units are assumed well mixed, so this sets the vertical resolution of grain-size composition that
    /// transport sees: it is a physical parameter. Missing means `record_dt_yr`. Must divide `record_dt_yr`;
    /// with it fixed, the record interval is a pure output choice.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layer_dt_yr: Option<f64>,
    /// Largest solid thickness (m) a unit may reach before a new one starts. Units are well mixed, so this
    /// bounds in thickness — independent of any clock — how coarsely transport sees composition with depth.
    /// Missing means [`DEFAULT_MAX_UNIT_SOLID_M`]; `null` means no cap. Always written out (as `null` when
    /// off), so a saved config means the same thing when read back.
    #[serde(default = "default_max_unit_solid")]
    pub max_unit_solid_m: Option<f64>,
    /// Flexural isostasy: the basement also moves under the changing weight of sediment and water. Missing
    /// means none (the basement moves only by `subsidence`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub isostasy: Option<Isostasy>,
    pub initial: Initial,
    pub sea_level: SeaLevel,
    pub subsidence: Subsidence,
    pub supply: Supply,
    pub transport: Transport,
    pub classes: Vec<GrainClass>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Initial {
    pub basement: Basement,
    /// Solid thickness of class 0 already present at start, everywhere. Zero means bare basement.
    #[serde(default)]
    pub sediment_solid_m: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Basement {
    /// Linear ramp from `top_m` at x=0 to `bottom_m` at x=length.
    Ramp { top_m: f64, bottom_m: f64 },
    /// One elevation per cell center.
    Profile(Vec<f64>),
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SeaLevel {
    pub base_m: f64,
    pub amplitude_m: f64,
    pub period_yr: f64,
    pub trend_m_per_yr: f64,
}

impl SeaLevel {
    pub fn at(&self, t_yr: f64) -> f64 {
        let osc = if self.period_yr > 0.0 {
            self.amplitude_m * (2.0 * std::f64::consts::PI * t_yr / self.period_yr).sin()
        } else {
            0.0
        };
        self.base_m + osc + self.trend_m_per_yr * t_yr
    }
}

/// Tectonic subsidence rate, linear from the proximal to the distal end. Positive lowers the basement.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Subsidence {
    pub proximal_m_per_yr: f64,
    pub distal_m_per_yr: f64,
}

impl Subsidence {
    pub fn rate(&self, x_frac: f64) -> f64 {
        self.proximal_m_per_yr + (self.distal_m_per_yr - self.proximal_m_per_yr) * x_frac
    }
}

/// Sediment entering through the proximal (x=0) boundary.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Supply {
    pub flux_m2_per_yr: f64,
    /// Fraction of the supply in each grain class. Must sum to 1.
    pub fractions: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Boundary {
    /// No sediment leaves the distal end.
    Closed,
    /// Sediment leaves the distal end at the slope of the last face (zero-gradient ghost cell would give no flux,
    /// so the last cell's own slope is reused).
    Open,
}

/// Diffusivity: `k_land` above sea level; below, `k_marine * exp(-depth / depth_decay_m) + k_deep`.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Transport {
    pub k_land_m2_per_yr: f64,
    pub k_marine_m2_per_yr: f64,
    pub depth_decay_m: f64,
    pub k_deep_m2_per_yr: f64,
    /// Solid thickness at the top of a column whose composition sets what leaves it.
    pub active_layer_m: f64,
    pub distal_boundary: Boundary,
    /// Fraction of the explicit stability limit used for the substep.
    pub cfl: f64,
    /// Width (m of water depth) of the land–marine diffusivity blend at the shoreline.
    /// 0 is a sharp switch at sea level; missing from older configs means 0.
    #[serde(default)]
    pub shoreline_width_m: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GrainClass {
    pub name: String,
    pub land_k_mult: f64,
    pub marine_k_mult: f64,
    /// Athy surface porosity.
    pub phi0: f64,
    /// Athy e-folding depth (m).
    pub compaction_length_m: f64,
}

impl Config {
    pub fn dx(&self) -> f64 {
        self.length_m / self.nx as f64
    }

    pub fn layer_dt(&self) -> f64 {
        self.layer_dt_yr.unwrap_or(self.record_dt_yr)
    }

    pub fn n_records(&self) -> usize {
        (self.duration_yr / self.record_dt_yr).round() as usize
    }

    /// Check the config before any computation.
    pub fn validate(&self) -> crate::Result<()> {
        self.check().map_err(crate::Error::Config)
    }

    // `!(x > 0.0)` is deliberate here: unlike `x <= 0.0` it is also true for NaN, so a NaN from the
    // config is rejected rather than passed through.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn check(&self) -> Result<(), String> {
        let nc = self.classes.len();
        if self.nx < 3 {
            return Err("nx must be at least 3".into());
        }
        if !(self.length_m > 0.0 && self.duration_yr > 0.0 && self.record_dt_yr > 0.0) {
            return Err("length, duration and record_dt must be positive".into());
        }
        if self.n_records() == 0 {
            return Err("duration must cover at least one record step".into());
        }
        // The run covers n_records × record_dt; anything else would silently run a different duration.
        let ratio = self.duration_yr / self.record_dt_yr;
        if (ratio - ratio.round()).abs() > 1e-9 * ratio.max(1.0) {
            return Err(format!(
                "duration_yr ({}) must be a whole multiple of record_dt_yr ({})",
                self.duration_yr, self.record_dt_yr
            ));
        }
        if nc == 0 {
            return Err("at least one grain class is required".into());
        }
        if self.supply.fractions.len() != nc {
            return Err(format!(
                "supply has {} fractions for {} classes",
                self.supply.fractions.len(),
                nc
            ));
        }
        let sum: f64 = self.supply.fractions.iter().sum();
        if (sum - 1.0).abs() > 1e-9 || self.supply.fractions.iter().any(|f| *f < 0.0) {
            return Err(format!(
                "supply fractions must be non-negative and sum to 1, got {sum}"
            ));
        }
        // Negative supply would drain the proximal cell past its stock and break conservation.
        if !(self.supply.flux_m2_per_yr >= 0.0) {
            return Err(format!(
                "supply.flux_m2_per_yr must be non-negative, got {}",
                self.supply.flux_m2_per_yr
            ));
        }
        // Names key per-class outputs (budget, PROV properties), where non-alphanumerics become '_'.
        let mut keys: Vec<String> = Vec::with_capacity(nc);
        for c in &self.classes {
            if !(0.0..1.0).contains(&c.phi0) || !(c.compaction_length_m > 0.0) {
                return Err(format!(
                    "class {}: phi0 must be in [0,1) and compaction length > 0",
                    c.name
                ));
            }
            // A negative multiplier is a negative diffusivity: transport runs uphill and is unstable.
            if !(c.land_k_mult >= 0.0) || !(c.marine_k_mult >= 0.0) {
                return Err(format!(
                    "class {}: land_k_mult and marine_k_mult must be non-negative",
                    c.name
                ));
            }
            let key: String = c
                .name
                .chars()
                .map(|ch| if ch.is_ascii_alphanumeric() { ch } else { '_' })
                .collect();
            if c.name.trim().is_empty() {
                return Err("class names must not be empty".into());
            }
            if keys.contains(&key) {
                return Err(format!(
                    "class names must be distinct once non-alphanumeric characters become '_': \"{}\" clashes",
                    c.name
                ));
            }
            keys.push(key);
        }
        if !(self.initial.sediment_solid_m >= 0.0) {
            return Err(format!(
                "initial.sediment_solid_m must be non-negative, got {}",
                self.initial.sediment_solid_m
            ));
        }
        if !(self.sea_level.period_yr >= 0.0) {
            return Err(format!(
                "sea_level.period_yr must be non-negative (0 = no oscillation), got {}",
                self.sea_level.period_yr
            ));
        }
        if let Basement::Profile(p) = &self.initial.basement {
            if p.len() != self.nx {
                return Err(format!(
                    "basement profile has {} values for nx={}",
                    p.len(),
                    self.nx
                ));
            }
        }
        let t = &self.transport;
        if t.k_land_m2_per_yr < 0.0 || t.k_marine_m2_per_yr < 0.0 || t.k_deep_m2_per_yr < 0.0 {
            return Err("diffusivities must be non-negative".into());
        }
        if !(t.shoreline_width_m >= 0.0) {
            return Err("shoreline_width_m must be non-negative".into());
        }
        if !(t.depth_decay_m > 0.0) || !(t.cfl > 0.0 && t.cfl <= 1.0) || !(t.active_layer_m > 0.0) {
            return Err("depth_decay_m and active_layer_m must be positive; cfl in (0,1]".into());
        }
        if let Some(l) = self.layer_dt_yr {
            let ratio = self.record_dt_yr / l;
            if !(l > 0.0) || (ratio - ratio.round()).abs() > 1e-9 || ratio.round() < 1.0 {
                return Err("layer_dt_yr must be positive and divide record_dt_yr".into());
            }
        }
        if !(self.resolution_m > 0.0) {
            return Err("resolution_m must be positive".into());
        }
        if let Some(iso) = &self.isostasy {
            iso.check()?;
        }
        // A cap of zero or less would start a new unit on every deposit, every substep.
        if let Some(cap) = self.max_unit_solid_m {
            if !(cap > 0.0) {
                return Err(format!(
                    "max_unit_solid_m must be positive when set, got {cap}"
                ));
            }
        }
        Ok(())
    }

    /// A prograding sand–mud wedge over one sea-level cycle and a slow rise. The reference configuration, used by the CLI.
    pub fn example() -> Self {
        Config {
            length_m: 100_000.0,
            nx: 200,
            duration_yr: 2_000_000.0,
            record_dt_yr: 20_000.0,
            start_age_ma: 12.0,
            resolution_m: 0.25,
            layer_dt_yr: None,
            max_unit_solid_m: default_max_unit_solid(),
            isostasy: None,
            initial: Initial {
                basement: Basement::Ramp {
                    top_m: 20.0,
                    bottom_m: -600.0,
                },
                sediment_solid_m: 0.0,
            },
            sea_level: SeaLevel {
                base_m: 0.0,
                amplitude_m: 30.0,
                period_yr: 800_000.0,
                trend_m_per_yr: 1.0e-5,
            },
            subsidence: Subsidence {
                proximal_m_per_yr: 2.0e-5,
                distal_m_per_yr: 1.2e-4,
            },
            supply: Supply {
                flux_m2_per_yr: 6.0,
                fractions: vec![0.4, 0.6],
            },
            transport: Transport {
                k_land_m2_per_yr: 3_000.0,
                k_marine_m2_per_yr: 800.0,
                depth_decay_m: 25.0,
                k_deep_m2_per_yr: 8.0,
                active_layer_m: 1.0,
                distal_boundary: Boundary::Closed,
                cfl: 0.4,
                shoreline_width_m: 5.0,
            },
            classes: vec![
                // Athy parameters after Sclater & Christie (1980).
                GrainClass {
                    name: "sand".into(),
                    land_k_mult: 1.0,
                    marine_k_mult: 0.35,
                    phi0: 0.49,
                    compaction_length_m: 3_700.0,
                },
                GrainClass {
                    name: "mud".into(),
                    land_k_mult: 1.0,
                    marine_k_mult: 1.6,
                    phi0: 0.63,
                    compaction_length_m: 1_960.0,
                },
            ],
        }
    }
}

/// Flexural isostasy of a thin elastic plate (Turcotte & Schubert, *Geodynamics*, ch. 3):
///
///   D w'''' + Δρ g w = q,   D = E Te³ / (12 (1 − ν²))
///
/// `w` is the downward deflection of the basement since the start, `q` the change since the start in the
/// weight of each column (grains, pore water, and the sea water above the bed, so sea-level change loads
/// the plate too). The deflection fills with sea water where the bed is under water (Δρ = ρ_mantle −
/// ρ_water) and with nothing where it is above (Δρ = ρ_mantle). Both ends of the section are free.
/// `subsidence` is read as water-loaded tectonic subsidence (the backstripping convention, Steckler & Watts
/// 1978): it already includes the response to its own water fill, so that water is not loaded again, and
/// where the bed is above sea level the missing water is an unloading (the column rebounds by
/// ρ_water T / ρ_mantle, the air-loaded equivalent of T). This keeps the equilibrium continuous where the
/// bed crosses sea level.
/// `elastic_thickness_m` 0 is local (Airy) isostasy. Densities default to 3300 (mantle), 1025 (sea water)
/// and 2650 kg/m³ (grains, all classes).
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Isostasy {
    pub elastic_thickness_m: f64,
    #[serde(default = "default_mantle_density")]
    pub mantle_density_kg_m3: f64,
    #[serde(default = "default_water_density")]
    pub water_density_kg_m3: f64,
    #[serde(default = "default_grain_density")]
    pub grain_density_kg_m3: f64,
    #[serde(default = "default_youngs_modulus")]
    pub youngs_modulus_pa: f64,
    #[serde(default = "default_poisson_ratio")]
    pub poisson_ratio: f64,
    /// E-folding time (yr) for the plate to reach equilibrium with its load (mantle flow). Missing means it
    /// responds at once.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub response_time_yr: Option<f64>,
}

/// Unit cap used when a config does not set `max_unit_solid_m`. On the reference configuration over 400 kyr, 0.5 m halves how
/// much the result depends on the layer clock (RMS surface difference between 20 kyr and 2 kyr clocks
/// 3.38 to 1.50 m) as well as 0.1 m does (1.51 m), at 1.75x the uncapped runtime against 7x.
pub const DEFAULT_MAX_UNIT_SOLID_M: f64 = 0.5;

fn default_max_unit_solid() -> Option<f64> {
    Some(DEFAULT_MAX_UNIT_SOLID_M)
}

fn default_mantle_density() -> f64 {
    3300.0
}
fn default_water_density() -> f64 {
    1025.0
}
fn default_grain_density() -> f64 {
    2650.0
}
fn default_youngs_modulus() -> f64 {
    7.0e10
}
fn default_poisson_ratio() -> f64 {
    0.25
}

impl Isostasy {
    /// Flexural rigidity D (N·m).
    pub fn rigidity(&self) -> f64 {
        let te = self.elastic_thickness_m;
        self.youngs_modulus_pa * te * te * te
            / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }

    // `!(x > 0.0)` also rejects NaN (see `Config::check`).
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn check(&self) -> Result<(), String> {
        if !(self.elastic_thickness_m >= 0.0) {
            return Err("isostasy.elastic_thickness_m must be non-negative".into());
        }
        if !(self.youngs_modulus_pa > 0.0) || !(0.0..0.5).contains(&self.poisson_ratio) {
            return Err(
                "isostasy: youngs_modulus_pa must be positive and poisson_ratio in [0, 0.5)".into(),
            );
        }
        let (m, w, g) = (
            self.mantle_density_kg_m3,
            self.water_density_kg_m3,
            self.grain_density_kg_m3,
        );
        if !(w > 0.0 && g > w && m > w) {
            return Err(format!(
                "isostasy densities must satisfy 0 < water < grain and water < mantle, got water {w}, grain {g}, mantle {m}"
            ));
        }
        if let Some(t) = self.response_time_yr {
            if !(t > 0.0) {
                return Err(format!(
                    "isostasy.response_time_yr must be positive when set, got {t}"
                ));
            }
        }
        Ok(())
    }
}
