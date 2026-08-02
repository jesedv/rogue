//! Exact breather solutions of the **standard-convention** focusing cubic NLS
//!
//! ```text
//! i u_t + u_xx + 2|u|² u = 0          (unit background u = e^{2it})
//! ```
//!
//! The unified breather formula (Akhmediev / Kuznetsov–Ma / Peregrine,
//! see Onorato et al., Phys. Rep. 2013):
//!
//! ```text
//! u(x,t) = [cos(2φ) cosh(Ωt) − i sin(2φ) sinh(Ωt) − cos(p x)]
//!        / [cosh(Ωt) − cos(φ) cos(p x)] · e^{2it}
//! ```
//! with `p = 2 sin(φ)`, `Ω = 2 sin(2φ)`. Real `φ` ⇒ **Akhmediev breather**
//! (periodic in x, breathes in t); imaginary `φ = iφ₀` ⇒ **Kuznetsov–Ma**
//! (localized in x, periodic in t); `φ → 0` limit ⇒ **Peregrine soliton**.
//!
//! These are solutions of the `u_xx + 2|u|²u` convention. Our solver uses
//! the convention `i ψ_t + ψ_xx + |ψ|²ψ = 0` (β2 = 2, focusing), which is
//! reached by the amplitude scaling `ψ(x,t) = √2·u(x,t)` (no time rescale;
//! `i·√2u_t + √2u_xx + |√2u|²·√2u = √2(iu_t + u_xx + 2|u|²u) = 0`). The
//! [`BreatherType`] mapping is exposed so tests and examples translate
//! exactly.

use num_complex::Complex64;
use rogue_nls::field::{Field1D, C};

/// Breather family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BreatherKind {
    Akhmediev,
    KuznetsovMa,
    Peregrine,
}

/// Type of a complex `phi` parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BreatherType {
    Akhmediev,
    KuznetsovMa,
    Peregrine,
}

/// Unified breather on the unit background `e^{2it}` (standard convention):
///
/// ```text
/// u(x,t) = [ cosh(Ωt − 2iφ) − cos(φ) cos(px) ] / [ cosh(Ωt) − cos(φ) cos(px) ] · e^{2it}
/// ```
/// with `p = 2 sin(φ)`, `Ω = 2 sin(2φ)`. Real `φ` ⇒ Akhmediev, `φ = iφ₀`
/// ⇒ Kuznetsov–Ma. `t_offset` shifts the breather's maximal-growth instant
/// to `t = -t_offset`.
pub fn breather(x: f64, t: f64, phi: Complex64, t_offset: f64) -> C {
    let p = 2.0 * phi.sin();
    let om = 2.0 * (2.0 * phi).sin();
    let tc = t + t_offset;
    let ch = (om * tc).cosh();
    let sh = (om * tc).sinh();
    let i = Complex64::new(0.0, 1.0);
    let cos_phi = phi.cos();
    let cos_px = (p * x).cos();
    // cosh(Ωt − 2iφ) = cos(2φ)cosh(Ωt) − i sin(2φ)sinh(Ωt)
    let num = (2.0 * phi).cos() * ch - i * (2.0 * phi).sin() * sh - cos_phi * cos_px;
    let den = ch - cos_phi * cos_px;
    let carrier = C::from_polar(1.0, 2.0 * t);
    num * carrier / den
}

/// Akhmediev breather, `phi` real in `(0, π/2)`. Perturbation period
/// `2π/p = π/sin φ`. Maximal growth at `t = -t_offset`.
pub fn akhmediev(x: f64, t: f64, phi: f64, t_offset: f64) -> C {
    breather(x, t, C::new(phi, 0.0), t_offset)
}

/// Kuznetsov–Ma breather, `phi0 > 0` (imaginary `φ = i φ₀`).
pub fn km_breather(x: f64, t: f64, phi0: f64) -> C {
    breather(x, t, C::new(0.0, phi0), 0.0)
}

/// Peregrine soliton: `u = e^{2it}(1 − 4(1+4it)/(1+4x²+16t²))`.
pub fn peregrine(x: f64, t: f64) -> C {
    let c = 1.0 + 4.0 * x * x + 16.0 * t * t;
    let num = C::new(4.0, 16.0 * t);
    let carrier = C::from_polar(1.0, 2.0 * t);
    (C::new(1.0, 0.0) - num.unscale(c)) * carrier
}

pub fn breather_type(phi: Complex64) -> BreatherType {
    if phi.re.abs() < 1e-9 && phi.im.abs() < 1e-9 {
        BreatherType::Peregrine
    } else if phi.im.abs() < 1e-12 {
        BreatherType::Akhmediev
    } else {
        BreatherType::KuznetsovMa
    }
}

/// Classification helper for the parameter given to the field constructors.
pub fn classify(phi: f64, is_km: bool) -> BreatherType {
    if phi.abs() < 1e-9 && !is_km {
        BreatherType::Peregrine
    } else if is_km {
        BreatherType::KuznetsovMa
    } else {
        BreatherType::Akhmediev
    }
}

/// Akhmediev initial condition for the solver at solver-time `t`.
///
/// `ψ(x,t) = √2·u(x,t)` with `u` the standard breather (solver convention
/// `iψ_t + ψ_xx + |ψ|²ψ = 0`). Box `lx` should be a multiple of the
/// perturbation period `π/sin φ`.
pub fn akhmediev_field(nx: usize, lx: f64, phi: f64, t: f64) -> Field1D {
    Field1D {
        data: (0..nx)
            .map(|i| {
                let x = i as f64 * lx / nx as f64;
                std::f64::consts::SQRT_2 * akhmediev(x, t, phi, 0.0)
            })
            .collect(),
        nx,
        lx,
    }
}

/// Peregrine initial condition for the solver at solver-time `t`.
pub fn peregrine_field(nx: usize, lx: f64, t: f64) -> Field1D {
    Field1D {
        data: (0..nx)
            .map(|i| {
                let x = i as f64 * lx / nx as f64;
                std::f64::consts::SQRT_2 * peregrine(x - lx / 2.0, t)
            })
            .collect(),
        nx,
        lx,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn peregrine_peak_three_times_background() {
        let c = peregrine(0.0, 0.0);
        assert!((c.norm() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn peregrine_decays_to_background() {
        let c = peregrine(50.0, 0.0);
        assert!((c.norm() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn akhmediev_is_periodic_in_x() {
        let phi: f64 = 0.8;
        let period = PI / phi.sin();
        let c0 = akhmediev(0.0, 1.0, phi, 0.0);
        let c1 = akhmediev(period, 1.0, phi, 0.0);
        assert!((c0 - c1).norm() < 1e-9);
    }

    #[test]
    fn akhmediev_amplitude_at_crest() {
        // At t=0, x=0 the AB amplitude is |1 + (cos2φ−1)/(1−cosφ)|.
        let phi = 1.0;
        let c = akhmediev(0.0, 0.0, phi, 0.0);
        let expected = (1.0 + ((2.0 * phi).cos() - 1.0) / (1.0 - phi.cos())).abs();
        assert!((c.norm() - expected).abs() < 1e-9);
    }

    #[test]
    fn km_breather_localized() {
        let c0 = km_breather(0.0, 0.0, 0.5);
        let c_far = km_breather(20.0, 0.0, 0.5);
        assert!(c0.norm() > 3.0);
        assert!((c_far.norm() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn classify_works() {
        assert_eq!(breather_type(C::new(1.0, 0.0)), BreatherType::Akhmediev);
        assert_eq!(breather_type(C::new(0.0, 1.0)), BreatherType::KuznetsovMa);
        assert_eq!(classify(0.0, false), BreatherType::Peregrine);
    }
}
