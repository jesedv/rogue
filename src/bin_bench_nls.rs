//! NLS performance benchmark. Verifies the hard constraint:
//! 2D step at 1024×1024 < 100 ms (release build).

use rogue_nls::params::NlsParams;
use rogue_nls::solver::{NlsSolver1D, NlsSolver2D};
use std::time::Instant;

fn bench_1d(nx: usize, steps: usize) {
    let params = NlsParams::focusing(20.0, 1e-3, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = rogue_nls::scenario::gaussian_1d(nx, 20.0, 1.0, 1.0, 10.0);
    // warmup
    for _ in 0..8 {
        solver.step(&mut field);
    }
    let t0 = Instant::now();
    for _ in 0..steps {
        solver.step(&mut field);
    }
    let dt = t0.elapsed().as_secs_f64() / steps as f64;
    println!(
        "1D nx={nx:<7} {steps:>5} steps  {:9.1} µs/step  ({:.1} Mpt/s)",
        dt * 1e6,
        nx as f64 / dt / 1e6
    );
}

fn bench_2d(nx: usize, ny: usize, steps: usize) {
    let params = NlsParams::focusing(20.0, 1e-3, 3.0);
    let mut solver = NlsSolver2D::new(params, nx, ny);
    let mut field = rogue_nls::scenario::gaussian_2d(nx, ny, 20.0, 20.0, 1.0, 2.0, 2.0, 10.0, 10.0);
    for _ in 0..2 {
        solver.step(&mut field);
    }
    let t0 = Instant::now();
    for _ in 0..steps {
        solver.step(&mut field);
    }
    let dt = t0.elapsed().as_secs_f64() / steps as f64;
    println!(
        "2D {nx}x{ny}  {steps:>3} steps  {:9.1} ms/step  ({:.1} Mpt/s)",
        dt * 1e3,
        (nx * ny) as f64 / dt / 1e6
    );
}

fn main() {
    println!("rogue-nls benchmark (release, threads={})", rayon::current_num_threads());
    bench_1d(4096, 1000);
    bench_1d(16384, 1000);
    bench_2d(256, 256, 50);
    bench_2d(512, 512, 20);
    bench_2d(1024, 1024, 10);
}
