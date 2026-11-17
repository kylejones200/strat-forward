//! `strat-forward`: a 2D forward stratigraphic model.
//!
//! - Accommodation: sea level (oscillation + trend) and hinged tectonic subsidence.
//! - Sediment mass balance: conservative finite volume, exact per-class solid budget.
//! - Transport: depth-dependent diffusion, separate diffusivity per grain class.
//! - Compaction: Athy porosity–depth law, recomputed at every record step.
//!
//! The record keeps deposition events (step, solids per class, water depth), not just the final grid,
//! so Wheeler diagrams, synthetic wells and cross sections come from the same run.

pub mod compaction;
pub mod config;
pub mod error;
pub mod export;
pub mod isostasy;
pub mod lineage;
pub mod model;
pub mod prov;
pub mod record;
pub mod reservoir;
pub mod rundir;
pub mod sequence;

pub use config::Config;
pub use error::{Error, Result};
pub use model::{History, Model};
