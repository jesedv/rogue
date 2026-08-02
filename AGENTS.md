# Rogue — Finite-Time Blow-up & Rogue Wave Predictor

## One-liner
Predict extreme wave events (rogue ocean waves, clear-air turbulence, plasma instabilities, fiber-optic damage) using the breakthrough mathematics of finite-time blow-up in supercritical defocusing NLS (Merle–Raphael–Rodnianski–Szeftel, 2026). For maritime, aerospace, plasma, and optics engineers.

## The Hard Math
- **Nonlinear Schrödinger equation (NLS)** — $i\partial_t \psi + \Delta \psi - |\psi|^{p-1}\psi = 0$.
- **Defocusing NLS** — supercritical case (Merle, Raphael, Rodnianski, Szeftel 2026).
- **Finite-time blow-up** — proven even in defocusing case.
- **Mass / energy / momentum criticality** — the subtle $L^2$-critical threshold.
- **Virial identities** — classical tool for blow-up.
- **Concentration-compactness** (P.-L. Lions, profile decomposition).
- **Soliton dynamics, modulation theory** — Grillakis, Weinstein, Schlag.
- **Compressible Euler / Navier–Stokes** — connection to fluid dynamics.
- **Rogue wave theory** — Akhmediev breathers, Kuznetsov–Ma breathers, Peregrine soliton.
- **Benjamin–Feir instability** — the basis of rogue wave formation.
- **Whitham modulation, inverse scattering**.

## The Real Problem
Rogue waves (extreme ocean waves, ≥ 2.2× significant wave height) cause:
- $2B+/yr in shipping damage.
- ~100 deaths/yr on cruise ships.
- Offshore platform damage.
- Submarine cable damage (communications).
- Wind turbine damage (offshore wind).

Current prediction: heuristic, statistical, no physics-based real-time.
- **WAVEWATCH III, SWAN** — spectral, no individual waves.
- **CFD (OpenFOAM, etc.)** — accurate but slow, no real-time.
- **Statistical (Rayleigh, Longuet-Higgins)** — underestimates extremes.

We use the **Merle–Raphael–Rodnianski–Szeftel (2026) breakthrough** to give:
- Mathematically rigorous blow-up detection.
- Real-time NLS solver for prediction.
- Connection to compressible Euler for fluid equivalence.
- Generalization to plasma, fiber optics, supercontinuum.

Other applications:
- **Clear-air turbulence** (aviation safety).
- **Plasma instabilities** (fusion reactor design, ionospheric phenomena).
- **Fiber-optic damage** (high-power laser, supercontinuum generation).
- **Bose–Einstein condensates** (BEC dynamics).

## Tech Stack
- **Rust** — NLS solver core.
- **`faer` / `ndarray`** — linear algebra.
- **WASM** — browser-based prediction.
- **CUDA / wgpu** — GPU acceleration for real-time.
- **TypeScript** — UI; **uPlot / Plotly / VTK.js** — wave-field viz.
- **Python bindings (PyO3)** — research adoption.
- **C** — SIMD inner loops.
- **Symplectic / geometric integrators** — for long-time behavior.

## Repository Layout
```
rogue/
├── Cargo.toml
├── crates/
│   ├── rogue-nls/          # NLS solver
│   ├── rogue-blow-up/      # blow-up detection
│   ├── rogue-rogue/        # rogue-wave detection (Akhmediev etc.)
│   ├── rogue-fluid/        # compressible Euler bridge
│   ├── rogue-plasma/       # plasma instabilities
│   ├── rogue-fiber/        # fiber-optic / supercontinuum
│   ├── rogue-turbulence/   # clear-air turbulence
│   ├── rogue-wasm/         # wasm-bindgen
│   └── rogue-production/  # real-data intake + physical forecast bridge
├── education/              # interactive browser demo (Vite + TS + WASM)
├── data/                   # sample sea-state CSVs for `rogue predict`
├── src/                    # CLI (akhmediev/peregrine/ocean/blowup/predict)
├── examples/               # Akhmediev, Peregrine, ocean data
└── docs/
    └── math.md
```

## Build & Test
- `cargo test`
- `cd education && npm run dev`
- `cargo bench`  (NLS solve time)
- `./dev p data/sample-buoy.csv`  (production forecast)
- `./scripts/regress-production.sh`  (production CSV intake regression)
- `./scripts/regress-akhmediev.sh`  (replay Akhmediev breathers)
- `./scripts/regress-peregrine.sh`  (replay Peregrine soliton)
- `./scripts/regress-ocean.sh`  (vs historical rogue wave events)
- `./scripts/regress-blowup.sh`  (verify blow-up detection on supercritical NLS)

## Conventions
- All algorithms reference the original paper.
- Energy / mass / momentum conservation reported.
- Reproducible (seeded).
- Bit-exact with published examples.
- Documented $L^2$ / $H^1$ regularity.

## Hard Constraints
- 2D NLS in < 100 ms for 1024×1024 grid.
- Blow-up detection: < 1 s for any $L^2$ norm.
- WASM bundle ≤ 8 MB.
- Numerical stability: symplectic / geometric integrators for long-time behavior.
- Reproducible (seeded).

## Non-Goals
- CFD for general flows (separate product, see atmosolver).
- Generic ocean modeling (WAVEWATCH III is good for spectral).
- ML-only prediction (v2: hybrid).
- Real-time vessel-mounted computation (edge devices).

## Open Questions
- Free vs. paid: free for personal, paid for shipping / offshore.
- White-label for marine insurers (Lloyd's, Allianz).
- Integration with existing ocean buoys / satellite data.
- Patent: novel blow-up detection algorithm.
- B2B vs. consumer: should we ship to vessel operators or insurers?

## References
- Merle, Raphael, Rodnianski, Szeftel, "Defocusing Wave Singularities" (2026, 2026 Breakthrough Prize in Mathematics).
- Akhmediev, Ankiewicz, *Nonlinear Pulses and Beams* (Chapman & Hall 1997).
- Akhmediev, Ankiewicz, Taki, "Waves that Appear from Nowhere and Disappear without a Trace" (Phys. Lett. A 2009) — Peregrine soliton.
- Onorato, Residori, Bortolozzo, et al., "Rogue Waves and Their Generating Mechanisms" (Phys. Rep. 2013).
- Dudley, Genty, Mussot, et al., "Rogue Waves and Analogies in Optics and Oceanography" (Nature Rev. Phys. 2019).
- Tao, *Nonlinear Dispersive Equations* (CBMS 2006).
- Weinstein, "Nonlinear Schrödinger Equations and Sharp Interpolation Estimates" (Comm. Math. Phys. 1983).
- Grillakis, "Regularity and Asymptotic Behavior of the Schrödinger Equation" (1988).
- Bourgain, *Global Solutions of Nonlinear Schrödinger Equations* (AMS 1999).
- P.-L. Lions, "The Concentration-Compactness Principle in the Calculus of Variations" (Ann. IHP 1984).
- Benjamin, Feir, "The Disintegration of Wave Trains on Deep Water" (J. Fluid Mech. 1967).
- Zakharov, "Stability of Periodic Waves of Finite Amplitude on the Surface of a Deep Fluid" (J. Appl. Mech. Tech. Phys. 1968).
- Hairer, Lubich, Wanner, *Geometric Numerical Integration* (Springer 2006).
- Peregrine, "Water Waves, Nonlinear Schrödinger Equations and Their Solutions" (J. Austral. Math. Soc. B 1983).
