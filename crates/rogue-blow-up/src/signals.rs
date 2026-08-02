use serde::{Deserialize, Serialize};

/// Instantaneous blow-up-relevant signals measured from a 1D/2D field.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct BlowUpSignals {
    /// Simulation time.
    pub t: f64,
    /// Conserved mass `∫|ψ|²`.
    pub mass: f64,
    /// Hamiltonian energy.
    pub energy: f64,
    /// `H¹` seminorm `‖∇ψ‖₂`.
    pub h1: f64,
    /// `‖∇ψ‖₂²` (the virial-relevant growth quantity).
    pub h1_sq: f64,
    /// Sup-norm `‖ψ‖_∞`.
    pub sup: f64,
    /// Virial about the center of mass `∫|x−x_c|²|ψ|²`.
    pub virial: f64,
    /// Center of mass (1D) or norm of 2D center offset.
    pub center: f64,
    /// True when the problem is `L²`-supercritical (p > 1 + 4/d).
    pub supercritical: bool,
    /// True when focusing and `E < 0` ⇒ rigorous finite-time blow-up.
    pub energy_negative: bool,
}

impl BlowUpSignals {
    /// Build from an NLS diagnostics snapshot (1D).
    pub fn from_diag_1d(
        t: f64,
        d: &rogue_nls::spectrum::Diagnostics,
        virial: f64,
        center: f64,
        power: f64,
        d_dim: usize,
    ) -> Self {
        Self {
            t,
            mass: d.mass,
            energy: d.energy,
            h1: d.h1,
            h1_sq: d.h1 * d.h1,
            sup: d.sup,
            virial,
            center,
            supercritical: power > 1.0 + 4.0 / d_dim as f64,
            energy_negative: d.energy < 0.0,
        }
    }
}

/// Ring buffer of recent blow-up signals (for growth-rate / ETA fitting).
pub struct SignalBuffer {
    capacity: usize,
    samples: Vec<BlowUpSignals>,
    head: usize,
    len: usize,
}

impl SignalBuffer {
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(2),
            samples: Vec::with_capacity(capacity.max(2)),
            head: 0,
            len: 0,
        }
    }

    pub fn push(&mut self, s: BlowUpSignals) {
        if self.len < self.capacity {
            self.samples.push(s);
            self.len += 1;
        } else {
            self.samples[self.head] = s;
            self.head = (self.head + 1) % self.capacity;
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Samples in chronological order.
    pub fn ordered(&self) -> Vec<BlowUpSignals> {
        (0..self.len)
            .map(|k| self.samples[(self.head + k) % self.capacity])
            .collect()
    }
}

/// Least-squares fit `y ≈ a + b x`, returns `(a, b, r²)`.
pub fn least_squares(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
    let n = xs.len() as f64;
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for (x, y) in xs.iter().zip(ys) {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-300 {
        return (0.0, 0.0, 0.0);
    }
    let b = (n * sxy - sx * sy) / denom;
    let a = (sy - b * sx) / n;
    let ybar = sy / n;
    let mut sse = 0.0;
    let mut sst = 0.0;
    for (x, y) in xs.iter().zip(ys) {
        let e = y - (a + b * x);
        sse += e * e;
        sst += (y - ybar) * (y - ybar);
    }
    let r2 = if sst > 0.0 { 1.0 - sse / sst } else { 0.0 };
    (a, b, r2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lsq_recovers_line() {
        let xs: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|x| 2.0 + 3.0 * x).collect();
        let (a, b, r2) = least_squares(&xs, &ys);
        assert!((a - 2.0).abs() < 1e-9);
        assert!((b - 3.0).abs() < 1e-9);
        assert!(r2 > 0.999999);
    }

    #[test]
    fn buffer_wraps() {
        let mut b = SignalBuffer::new(4);
        for i in 0..10 {
            let mut s = BlowUpSignals::default();
            s.t = i as f64;
            b.push(s);
        }
        assert_eq!(b.len(), 4);
        let ordered = b.ordered();
        let ts: Vec<f64> = ordered.iter().map(|s| s.t).collect();
        assert_eq!(ts, vec![6.0, 7.0, 8.0, 9.0]);
    }
}
