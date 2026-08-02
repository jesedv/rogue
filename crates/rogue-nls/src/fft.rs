use num_complex::Complex64;
use rayon::prelude::*;
use rustfft::{Fft, FftPlanner};
use std::cell::RefCell;
use std::sync::Arc;

pub type C = Complex64;

thread_local! {
    static SCRATCH: RefCell<Vec<C>> = const { RefCell::new(Vec::new()) };
}

/// 1D FFT engine with precomputed physical wavenumbers.
pub struct FftEngine1D {
    plan: Arc<dyn Fft<f64>>,
    inv_plan: Arc<dyn Fft<f64>>,
    /// Physical wavenumber `k_j = 2π·m_j/lx` with `m_j` in `0..nx/2` for
    /// the first half and `(j−nx)` for the second (DFT bin order).
    pub k: Vec<f64>,
    /// `k_j²`.
    pub ksq: Vec<f64>,
    lx: f64,
}

impl FftEngine1D {
    pub fn new(nx: usize, lx: f64) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        let plan = planner.plan_fft_forward(nx);
        let inv_plan = planner.plan_fft_inverse(nx);
        let k: Vec<f64> = (0..nx)
            .map(|j| {
                let m = if j <= nx / 2 { j as i64 } else { j as i64 - nx as i64 };
                2.0 * std::f64::consts::PI * m as f64 / lx
            })
            .collect();
        let ksq = k.iter().map(|v| v * v).collect();
        Self { plan, inv_plan, k, ksq, lx }
    }

#[inline]
pub fn forward(&self, data: &mut [C]) {
    let plan = &self.plan;
    SCRATCH.with(|s| {
        let mut sc = s.borrow_mut();
        sc.resize(self.plan.get_inplace_scratch_len(), C::default());
        plan.process_with_scratch(data, &mut sc);
    });
}

#[inline]
pub fn inverse(&self, data: &mut [C]) {
    let plan = &self.inv_plan;
    SCRATCH.with(|s| {
        let mut sc = s.borrow_mut();
        sc.resize(self.inv_plan.get_inplace_scratch_len(), C::default());
        plan.process_with_scratch(data, &mut sc);
    });
    // unnormalized rustfft → normalize by 1/n
    let inv = 1.0 / self.k.len() as f64;
    for c in data.iter_mut() {
        *c = c.scale(inv);
    }
}

    /// Linear half-step: `ψ̂ *= exp(−i (β2/2) k² (dt/2))`.
    pub fn apply_linear_half(&self, data: &mut [C], beta2: f64, dt: f64) {
        let phase_scale = -0.5 * beta2 * dt * 0.5;
        for (c, k2) in data.iter_mut().zip(&self.ksq) {
            let (s, co) = (phase_scale * k2).sin_cos();
            *c = C::new(c.re * co - c.im * s, c.re * s + c.im * co);
        }
    }

    /// `H¹` seminorm `‖∇ψ‖_{L²}` via Parseval:
    /// `sqrt( (lx/nx²) · Σ k² |ψ̂|² )`.
    pub fn h1_seminorm(&self, data: &[C]) -> f64 {
        let s2: f64 = data
            .iter()
            .zip(&self.ksq)
            .map(|(c, k2)| k2 * c.norm_sqr())
            .sum();
        let nx = self.k.len() as f64;
        (s2 * self.lx / (nx * nx)).sqrt()
    }
}

/// 2D FFT engine. Uses row/column passes with a transposition so every
/// transform works on contiguous memory; rows are processed in parallel.
pub struct FftEngine2D {
    row_plan: Arc<dyn Fft<f64>>,
    col_plan: Arc<dyn Fft<f64>>,
    inv_row_plan: Arc<dyn Fft<f64>>,
    inv_col_plan: Arc<dyn Fft<f64>>,
    pub nx: usize,
    pub ny: usize,
    pub kx: Vec<f64>,
    pub ky: Vec<f64>,
    /// `kx[i]² + ky[j]²` in row-major layout.
    pub ksq: Vec<f64>,
    lx: f64,
    ly: f64,
    tmp: Vec<C>,
}

impl FftEngine2D {
    pub fn new(nx: usize, ny: usize, lx: f64, ly: f64) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        let row_plan = planner.plan_fft_forward(nx);
        let col_plan = planner.plan_fft_forward(ny);
        let inv_row_plan = planner.plan_fft_inverse(nx);
        let inv_col_plan = planner.plan_fft_inverse(ny);
        let kx: Vec<f64> = (0..nx)
            .map(|j| {
                let m = if j <= nx / 2 { j as i64 } else { j as i64 - nx as i64 };
                2.0 * std::f64::consts::PI * m as f64 / lx
            })
            .collect();
        let ky: Vec<f64> = (0..ny)
            .map(|j| {
                let m = if j <= ny / 2 { j as i64 } else { j as i64 - ny as i64 };
                2.0 * std::f64::consts::PI * m as f64 / ly
            })
            .collect();
        let mut ksq = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                ksq[j * nx + i] = kx[i] * kx[i] + ky[j] * ky[j];
            }
        }
        let tmp = vec![C::default(); nx * ny];
        Self {
            row_plan,
            col_plan,
            inv_row_plan,
            inv_col_plan,
            nx,
            ny,
            kx,
            ky,
            ksq,
            lx,
            ly,
            tmp,
        }
    }

    fn rows(&self, data: &mut [C], plan: &Arc<dyn Fft<f64>>) {
        let slen = plan.get_inplace_scratch_len();
        let plan = plan.clone();
        data.par_chunks_mut(self.nx).for_each(|row| {
            SCRATCH.with(|s| {
                let mut sc = s.borrow_mut();
                sc.resize(slen, C::default());
                plan.process_with_scratch(row, &mut sc);
            });
        });
    }

    fn cols(&mut self, data: &mut [C], plan: &Arc<dyn Fft<f64>>) {
        let (nx, ny) = (self.nx, self.ny);
        let slen = plan.get_inplace_scratch_len();
        let plan = plan.clone();
        // transpose into tmp: tmp[i*ny + j] = data[j*nx + i]
        {
            let tmp = &mut self.tmp;
            tmp.par_chunks_mut(ny).enumerate().for_each(|(i, chunk)| {
                for j in 0..ny {
                    chunk[j] = data[j * nx + i];
                }
            });
        }
        // fft along former columns (now contiguous rows of length ny)
        {
            let tmp = &mut self.tmp;
            tmp.par_chunks_mut(ny).for_each(|row| {
                SCRATCH.with(|s| {
                    let mut sc = s.borrow_mut();
                    sc.resize(slen, C::default());
                    plan.process_with_scratch(row, &mut sc);
                });
            });
        }
        // transpose back
        {
            let tmp = &self.tmp;
            data.par_chunks_mut(nx).enumerate().for_each(|(j, row)| {
                for i in 0..nx {
                    row[i] = tmp[i * ny + j];
                }
            });
        }
    }

    pub fn forward(&mut self, data: &mut [C]) {
        let (row_plan, col_plan) = (self.row_plan.clone(), self.col_plan.clone());
        self.rows(data, &row_plan);
        self.cols(data, &col_plan);
    }

    pub fn inverse(&mut self, data: &mut [C]) {
        let (row_plan, col_plan) = (self.inv_row_plan.clone(), self.inv_col_plan.clone());
        self.rows(data, &row_plan);
        self.cols(data, &col_plan);
        let inv = 1.0 / (self.nx * self.ny) as f64;
        for c in data.iter_mut() {
            *c = c.scale(inv);
        }
    }

    /// Linear half-step: `ψ̂ *= exp(−i (β2/2) k² (dt/2))`.
    pub fn apply_linear_half(&self, data: &mut [C], beta2: f64, dt: f64) {
        let phase_scale = -0.5 * beta2 * dt * 0.5;
        for (c, k2) in data.iter_mut().zip(&self.ksq) {
            let (s, co) = (phase_scale * k2).sin_cos();
            *c = C::new(c.re * co - c.im * s, c.re * s + c.im * co);
        }
    }

    /// `H¹` seminorm via Parseval: `sqrt( (lx·ly)/(nx·ny)² · Σ k²|ψ̂|² )`.
    pub fn h1_seminorm(&self, data: &[C]) -> f64 {
        let s2: f64 = data
            .iter()
            .zip(&self.ksq)
            .map(|(c, k2)| k2 * c.norm_sqr())
            .sum();
        let n = (self.nx * self.ny) as f64;
        (s2 * self.lx * self.ly / (n * n)).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_1d() {
        let nx = 32;
        let mut eng = FftEngine1D::new(nx, 10.0);
        let mut buf: Vec<C> = (0..nx)
            .map(|i| C::new((i as f64 * 0.7).sin() + 0.3, (i as f64 * 0.4).cos()))
            .collect();
        let orig = buf.clone();
        eng.forward(&mut buf);
        eng.inverse(&mut buf);
        for (a, b) in buf.iter().zip(&orig) {
            assert!((a - b).norm() < 1e-9);
        }
    }

    #[test]
    fn roundtrip_2d() {
        let (nx, ny) = (16, 20);
        let mut eng = FftEngine2D::new(nx, ny, 5.0, 6.0);
        let mut buf: Vec<C> = (0..nx * ny)
            .map(|idx| {
                let i = idx % nx;
                let j = idx / nx;
                C::new(
                    (i as f64 * 0.3 + j as f64 * 0.2).sin(),
                    (i as f64 * 0.1 - j as f64 * 0.4).cos(),
                )
            })
            .collect();
        let orig = buf.clone();
        eng.forward(&mut buf);
        eng.inverse(&mut buf);
        for (a, b) in buf.iter().zip(&orig) {
            assert!((a - b).norm() < 1e-9);
        }
    }

    #[test]
    fn h1_of_plane_wave_is_exact() {
        // ψ = exp(i k0 x) with k0 on the discrete grid has |∇ψ|² = k0² M.
        let (nx, lx) = (64, 8.0);
        let k0 = 2.0 * std::f64::consts::PI * 3.0 / lx;
        let eng = FftEngine1D::new(nx, lx);
        let mut buf: Vec<C> = (0..nx)
            .map(|i| C::from_polar(1.0, k0 * i as f64 * lx / nx as f64))
            .collect();
        eng.forward(&mut buf);
        let h1 = eng.h1_seminorm(&buf);
        let expected = k0 * (nx as f64 * lx / nx as f64).sqrt();
        assert!((h1 - expected).abs() < 1e-9);
    }
}
