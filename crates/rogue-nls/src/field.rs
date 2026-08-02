use num_complex::Complex64;
use serde::{Deserialize, Serialize};

pub type C = Complex64;

/// 1D complex field on a periodic box `[0, lx)`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Field1D {
    pub data: Vec<C>,
    pub nx: usize,
    pub lx: f64,
}

impl Field1D {
    pub fn new(nx: usize, lx: f64) -> Self {
        Self {
            data: vec![C::new(0.0, 0.0); nx],
            nx,
            lx,
        }
    }

    #[inline]
    pub fn dx(&self) -> f64 {
        self.lx / self.nx as f64
    }

    #[inline]
    pub fn at(&self, i: usize) -> C {
        self.data[i]
    }

    /// |ψ|² at index `i`.
    #[inline]
    pub fn density(&self, i: usize) -> f64 {
        let c = self.data[i];
        c.norm_sqr()
    }

    pub fn sup_norm(&self) -> f64 {
        self.data.iter().map(|c| c.norm()).fold(0.0, f64::max)
    }

    pub fn mass(&self) -> f64 {
        self.data.iter().map(|c| c.norm_sqr()).sum::<f64>() * self.dx()
    }

    pub fn center_of_mass(&self) -> f64 {
        let m = self.mass();
        if m <= 0.0 {
            return 0.0;
        }
        let mut s = 0.0;
        for (i, c) in self.data.iter().enumerate() {
            let x = (i as f64 + 0.5) * self.dx();
            s += x * c.norm_sqr();
        }
        s * self.dx() / m
    }

    /// Virial `∫ (x − x_c)² |ψ|² dx` about the center of mass.
    pub fn virial_about_com(&self) -> f64 {
        let xc = self.center_of_mass();
        let mut s = 0.0;
        for (i, c) in self.data.iter().enumerate() {
            let x = (i as f64 + 0.5) * self.dx();
            s += (x - xc) * (x - xc) * c.norm_sqr();
        }
        s * self.dx()
    }

    /// Apply nonlinear phase `ψ *= exp(−i γ dt |ψ|^(p−1))` in place.
    #[inline]
    pub fn apply_nonlinear_phase(&mut self, gamma: f64, power: f64, dt: f64) {
        let exp = (power - 1.0) * 0.5;
        for c in self.data.iter_mut() {
            let a2 = c.norm_sqr();
            let phase = -gamma * a2.powf(exp) * dt;
            let (s, c_ph) = phase.sin_cos();
            *c = C::new(c.re * c_ph - c.im * s, c.re * s + c.im * c_ph);
        }
    }
}

/// 2D complex field on a periodic box `[0, lx) × [0, ly)` (row-major,
/// index `j * nx + i`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Field2D {
    pub data: Vec<C>,
    pub nx: usize,
    pub ny: usize,
    pub lx: f64,
    pub ly: f64,
}

impl Field2D {
    pub fn new(nx: usize, ny: usize, lx: f64, ly: f64) -> Self {
        Self {
            data: vec![C::new(0.0, 0.0); nx * ny],
            nx,
            ny,
            lx,
            ly,
        }
    }

    #[inline]
    pub fn dx(&self) -> f64 {
        self.lx / self.nx as f64
    }

    #[inline]
    pub fn dy(&self) -> f64 {
        self.ly / self.ny as f64
    }

    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    #[inline]
    pub fn at(&self, i: usize, j: usize) -> C {
        self.data[self.index(i, j)]
    }

    #[inline]
    pub fn density(&self, i: usize, j: usize) -> f64 {
        self.at(i, j).norm_sqr()
    }

    pub fn sup_norm(&self) -> f64 {
        self.data.iter().map(|c| c.norm()).fold(0.0, f64::max)
    }

    pub fn mass(&self) -> f64 {
        self.data.iter().map(|c| c.norm_sqr()).sum::<f64>() * self.dx() * self.dy()
    }

    pub fn center_of_mass(&self) -> (f64, f64) {
        let m = self.mass();
        if m <= 0.0 {
            return (0.0, 0.0);
        }
        let (mut sx, mut sy) = (0.0, 0.0);
        for j in 0..self.ny {
            let y = (j as f64 + 0.5) * self.dy();
            for i in 0..self.nx {
                let x = (i as f64 + 0.5) * self.dx();
                let d = self.at(i, j).norm_sqr();
                sx += x * d;
                sy += y * d;
            }
        }
        (sx * self.dx() * self.dy() / m, sy * self.dx() * self.dy() / m)
    }

    pub fn virial_about_com(&self) -> f64 {
        let (xc, yc) = self.center_of_mass();
        let mut s = 0.0;
        for j in 0..self.ny {
            let y = (j as f64 + 0.5) * self.dy();
            for i in 0..self.nx {
                let x = (i as f64 + 0.5) * self.dx();
                s += ((x - xc) * (x - xc) + (y - yc) * (y - yc)) * self.at(i, j).norm_sqr();
            }
        }
        s * self.dx() * self.dy()
    }

    #[inline]
    pub fn apply_nonlinear_phase(&mut self, gamma: f64, power: f64, dt: f64) {
        let exp = (power - 1.0) * 0.5;
        for c in self.data.iter_mut() {
            let a2 = c.norm_sqr();
            let phase = -gamma * a2.powf(exp) * dt;
            let (s, c_ph) = phase.sin_cos();
            *c = C::new(c.re * c_ph - c.im * s, c.re * s + c.im * c_ph);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonlinear_phase_preserves_mass() {
        let mut f = Field1D::new(64, 10.0);
        for (i, c) in f.data.iter_mut().enumerate() {
            let x = i as f64;
            *c = C::new((x * 0.3).sin() + 1.0, (x * 0.1).cos());
        }
        let m0 = f.mass();
        f.apply_nonlinear_phase(1.0, 3.0, 0.01);
        assert!((f.mass() - m0).abs() < 1e-12);
    }
}
