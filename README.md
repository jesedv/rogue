# Rogue — finite-time blow-up & rogue-wave predictor

Predict extreme wave events — rogue ocean waves, clear-air turbulence, plasma
instabilities, fiber-optic damage — from the mathematics of finite-time
blow-up in supercritical **defocusing** NLS
(Merle–Raphael–Rodnianski–Szeftel, 2026 Breakthrough Prize in Mathematics).

Real-time split-step NLS solver in Rust (→ WASM), with mathematically rigorous
blow-up detection, exact breather physics, a browser dashboard, **and a
production sea-state intake** that turns real observational data into physical
forecasts.

## Layout

```
crates/
  rogue-nls          # split-step NLS solver (1D/2D), FFT, diagnostics, scenarios
  rogue-blow-up      # virial/energy + H¹-acceleration blow-up detection, ETA
  rogue-rogue        # Akhmediev/KM/Peregrine breathers, Benjamin–Feir, crest detection
  rogue-fluid        # Madelung ↔ compressible-Euler bridge
  rogue-plasma       # plasma modulational instability
  rogue-fiber        # fiber-optic soliton power / bandwidth (supercontinuum)
  rogue-turbulence   # clear-air turbulence intermittency
  rogue-wasm         # wasm-bindgen bindings for the browser dashboard
  rogue-production   # real-data intake + physical forecast bridge (production)
education/           # interactive browser demo (Vite + TS + WASM)
data/                # sample sea-state CSV records for `rogue predict`
docs/math.md         # the full mathematical basis
scripts/             # regression + bench suite (incl. production regress)
```

## Quickstart

```bash
./dev setup     # cargo fetch + npm install + wasm-pack build
./dev test      # cargo test --workspace
./dev bench     # NLS performance (1024² must be < 100 ms/step)
./dev ui        # browser education demo at localhost:5173
./dev regress    # physics regression suite
./dev regress-prod # production intake regression
./dev r         # CLI help
```

## Production mode (real data → forecast)

`rogue-production` ingests real, dimensional sea-state observations and runs
the physical NLS forecast bridge. Input is a CSV/TSV with `t,hs,tp` (and
optional `gamma`):

```
t,hs,tp,gamma
0,1.5,8.0,3.3
3600,2.8,9.5,3.3
```

```bash
./dev p data/sample-buoy.csv        # human-readable report + per-obs forecast
./dev p data/sample-buoy.csv --json # machine-readable JSON report
```

The physical bridge ([`crates/rogue-production/src/sea.rs`]) derives the
carrier wavenumber from `Tp` via deep-water dispersion `ω² = g·k`, then the
group-velocity dispersion `β = ω₀/(8k₀²)` and nonlinearity `γ = ω₀k₀²/2` used
by the split-step solver. Output is dimensional: wave events have crest
factor + significant wave height in metres.

## Education / interactive

```bash
./dev wasm   # rebuild WASM into education/pkg
./dev ui     # run the browser playground (Akhmediev, Peregrine, JONSWAP ocean,
             # blow-up, soliton) rendered as a live 2D field
```

## CLI

```bash
cargo run --release -- akhmediev 1.0 1.0    # replay Akhmediev breather
cargo run --release -- peregrine 1.0        # Peregrine soliton evolution
cargo run --release -- ocean 60.0           # JONSWAP random sea
cargo run --release -- blowup 5.0 0.3       # quintic blow-up run + detector
cargo run --release -- soliton 4.0          # soliton walk
cargo run --release -- fiber 5.0            # fiber-optic diagnostics
cargo run --release -- predict data/sample-buoy.csv [--json]   # production
```

## Regression suite

```bash
./dev regress ./dev regress-prod
```

## Hard constraints (release build, 12 threads)

| Constraint | Target | Measured |
|---|---|---|
| 2D step 1024×1024 | < 100 ms | 59 ms |
| Blow-up detect | < 1 s | streaming fit |
| WASM bundle | ≤ 8 MB | verified at build |

## Math

`docs/math.md` — governing equations, invariants, Strang splitting, blow-up
detection (virial + H¹ ETA), Benjamin–Feir law `γ(ν)=ν√(2A²−ν²)`, unified
breather formula, Madelung/Euler bridge.

## License

Free and open for maritime, aviation, plasma, and optics use — this tool is
built to help any agency, company, or operator worldwide. MIT.