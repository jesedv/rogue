use crate::field::{Field1D, Field2D, C};
use crate::fft::{FftEngine1D, FftEngine2D};
use crate::params::NlsParams;
use serde::{Deserialize, Serialize};

/// Conserved quantities and regularity monitors, reported per the project
/// conventions (energy / mass / momentum conservation, L² and H¹ regularity).
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct Diagnostics {
    pub mass: f64,
    pub momentum: f64,
    pub energy: f64,
    pub h1: f64,
    pub sup: f64,
}

impl Diagnostics {
    /// Monitor on a 1D field. Performs one forward FFT internally.
    pub fn measure_1d(eng: &FftEngine1D, field: &Field1D, params: &NlsParams) -> Self {
        let mut spec = field.data.clone();
        eng.forward(&mut spec);
        Self::from_spectrum_1d(eng, &spec, &field.data, params)
    }

    /// Monitor on a field already given in spectral form (row-major 1D).
    pub fn from_spectrum_1d(
        eng: &FftEngine1D,
        spec: &[C],
        phys: &[C],
        params: &NlsParams,
    ) -> Self {
        let nx = eng.k.len();
        let dx = params.lx / nx as f64;
        let mass = phys.iter().map(|c| c.norm_sqr()).sum::<f64>() * dx;
        // P = (lx/nx²)·Σ k|ψ̂|²  (physical momentum ∫Im(conj ψ·ψ_x)dx)
        let momentum = spec
            .iter()
            .zip(&eng.k)
            .map(|(c, k)| k * c.norm_sqr())
            .sum::<f64>()
            * dx
            / nx as f64;
        let h1 = eng.h1_seminorm(spec);
        let pot: f64 = phys
            .iter()
            .map(|c| c.norm_sqr().powf((params.power + 1.0) * 0.5))
            .sum::<f64>()
            * dx;
        let energy = 0.25 * params.beta2 * h1 * h1 + params.gamma() * pot / (params.power + 1.0);
        let sup = phys.iter().map(|c| c.norm()).fold(0.0, f64::max);
        Self { mass, momentum, energy, h1, sup }
    }

    /// Monitor on a 2D field (row-major spectral buffer of length nx·ny).
    pub fn measure_2d(eng: &mut FftEngine2D, field: &Field2D, params: &NlsParams) -> Self {
        let mut spec = field.data.clone();
        eng.forward(&mut spec);
        Self::from_spectrum_2d(eng, &spec, &field.data, params)
    }

    pub fn from_spectrum_2d(
        eng: &FftEngine2D,
        spec: &[C],
        phys: &[C],
        params: &NlsParams,
    ) -> Self {
        let (nx, ny) = (eng.nx, eng.ny);
        let dx = params.lx / nx as f64;
        let dy = params.lx / ny as f64;
        let dv = dx * dy;
        let mass = phys.iter().map(|c| c.norm_sqr()).sum::<f64>() * dv;
        let momentum = spec
            .iter()
            .enumerate()
            .map(|(idx, c)| {
                let i = idx % nx;
                let j = idx / nx;
                (eng.kx[i] + eng.ky[j]) * c.norm_sqr()
            })
            .sum::<f64>()
            * dv
            / (nx * ny) as f64;
        let h1 = eng.h1_seminorm(spec);
        let pot: f64 = phys
            .iter()
            .map(|c| c.norm_sqr().powf((params.power + 1.0) * 0.5))
            .sum::<f64>()
            * dv;
        let energy = 0.25 * params.beta2 * h1 * h1 + params.gamma() * pot / (params.power + 1.0);
        let sup = phys.iter().map(|c| c.norm()).fold(0.0, f64::max);
        Self { mass, momentum, energy, h1, sup }
    }
}
