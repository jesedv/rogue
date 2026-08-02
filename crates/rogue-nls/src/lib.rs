//! Core nonlinear Schrödinger equation (NLS) solver.
//!
//! Equation (generalized, defocusing/focusing):
//! ```text
//! i ∂t ψ = -(β2/2) ∇²ψ + γ |ψ|^(p-1) ψ
//! ```
//! with dispersion `β2`, nonlinearity strength `γ` (`γ > 0` defocusing,
//! `γ < 0` focusing) and power `p`.
//!
//! The `defocusing`/`focusing` helpers build parameters matching the
//! Merle–Raphael–Rodnianski–Szeftel convention
//! `i ∂t ψ + Δψ − |ψ|^(p-1) ψ = 0` (i.e. `β2 = 2`, `γ = ±1`).
//!
//! The default integrator is a **second-order Strang split-step**
//! (linear half-step, nonlinear step, linear half-step), which is
//! symplectic and hence stable over long times.

pub mod fft;
pub mod field;
pub mod params;
pub mod random;
pub mod scenario;
pub mod solver;
pub mod spectrum;

pub use field::{Field1D, Field2D};
pub use params::{Focus, NlsParams};
pub use solver::{NlsSolver1D, NlsSolver2D};
