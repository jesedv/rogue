//! Compressible-Euler bridge for the NLS solver (Madelung transformation)
//! and a reference isothermal Euler solver for cross-validation.

pub mod euler;
pub mod madelung;

pub use euler::EulerState;
pub use madelung::{div_rho_u, madelung, MadelungState};
