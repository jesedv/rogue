use serde::{Deserialize, Serialize};

/// Concavity of the virial `V(t) = ∫|x−x_c|²|ψ|² dx`.
///
/// For focusing NLS, blow-up is classically signalled by `V″ < 0` (virial
/// collapsing), and the virial identity `V″ = 16E` (about the center of
/// mass, focusing cubic) links it directly to the energy criterion.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct VirialCurvature {
    pub d2: f64,
    pub sign: i8,
}

/// Finite-difference second derivative of the virial from a time series.
/// Uses central differences on the last three samples; returns `None` when
/// fewer than three samples are available.
pub fn virial_second_derivative(ts: &[f64], virial: &[f64]) -> Option<VirialCurvature> {
    if ts.len() < 3 {
        return None;
    }
    let n = ts.len();
    let (t0, t1, t2) = (ts[n - 3], ts[n - 2], ts[n - 1]);
    let (v0, v1, v2) = (virial[n - 3], virial[n - 2], virial[n - 1]);
    let h1 = t1 - t0;
    let h2 = t2 - t1;
    if h1.abs() < 1e-12 || h2.abs() < 1e-12 {
        return None;
    }
    // Non-uniform central second difference.
    let d2 = 2.0 * (v0 / (h1 * (h1 + h2)) - v1 / (h1 * h2) + v2 / (h2 * (h1 + h2)));
    Some(VirialCurvature {
        d2,
        sign: if d2 < 0.0 { -1 } else if d2 > 0.0 { 1 } else { 0 },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concave_down_virial_detected() {
        // V = -(t-2)^2 + 10 → concave down everywhere.
        let ts: Vec<f64> = (0..10).map(|i| 0.1 * i as f64).collect();
        let virial: Vec<f64> = ts.iter().map(|t| -((t - 2.0).powi(2)) + 10.0).collect();
        let c = virial_second_derivative(&ts, &virial).unwrap();
        assert_eq!(c.sign, -1);
        assert!(c.d2 < 0.0);
        assert!((c.d2 - (-2.0)).abs() < 0.5, "d2={}", c.d2);
    }

    #[test]
    fn linear_virial_flat() {
        let ts: Vec<f64> = (0..10).map(|i| 0.1 * i as f64).collect();
        let virial: Vec<f64> = ts.iter().map(|t| 3.0 * t + 1.0).collect();
        let c = virial_second_derivative(&ts, &virial).unwrap();
        assert!(c.d2.abs() < 1e-6);
    }
}
