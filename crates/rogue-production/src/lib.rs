//! rogue-production — real-world intake and physical forecast bridge.
//!
//! This crate turns real, dimensional sea-state data (Hs, Tp) into a physical
//! NLS forecast: it ingests CSV observations, derives the carrier-scaling
//! coefficients, runs the split-step solver, and emits a machine-readable
//! report of rogue events and finite-time-blow-up risk.
//!
//! Layout:
//! - [`sea`] — dimensional sea state and NLS scaling bridge.
//! - [`io`] — CSV / TSV intake of observations.
//! - [`forecast`] — the prediction pipeline + report.

pub mod forecast;
pub mod io;
pub mod sea;

pub use forecast::{forecast_observations, ForecastConfig, ForecastReport};
pub use io::{parse_observations, Observation};
pub use sea::SeaState;