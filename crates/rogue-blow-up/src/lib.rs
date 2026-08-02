//! Finite-time blow-up detection for supercritical focusing NLS.
//!
//! Rigorous inputs (Merle–Raphael–Rodnianski–Szeftel):
//! - negative energy `E < 0` with sub-threshold mass ⇒ guaranteed
//!   finite-time blow-up (virial argument);
//! - `L²`-supercritical regime `p > 1 + 4/d` ⇒ blow-up for negative energy.
//!
//! The runtime detector combines the energy criterion with growth-rate and
//! concentration signals and produces a least-squares estimate of the
//! blow-up time `T_*` from the self-similar law `‖∇ψ‖² ≈ C (T_* − t)⁻¹`.

mod detect;
mod signals;
mod virial;

pub use detect::{BlowUpDetector, BlowUpState};
pub use signals::{BlowUpSignals, SignalBuffer};
pub use virial::{virial_second_derivative, VirialCurvature};
