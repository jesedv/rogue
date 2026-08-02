<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="education/public/waves-white.svg">
    <img src="education/public/waves-green.svg" width="64" alt="Rogue" />
  </picture>
</p>

<h1 align="center">Rogue — finite-time blow-up & rogue-wave predictor</h1>

<p align="center">
  <strong>Real-time physics-based extreme-wave forecast.</strong><br />
  Ocean · aerospace · plasma · fiber-optics — NLS solver with mathematically rigorous blow-up detection.
</p>

<p align="center">
  <a href="https://github.com/jesedv/rogue/actions"><img src="https://img.shields.io/github/actions/workflow/status/jesed/rogue/ci.yml?branch=main" alt="CI" /></a>
  <a href="https://crates.io/crates/rogue-nls"><img src="https://img.shields.io/crates/v/rogue-nls" alt="crates.io" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="MIT" /></a>
  <a href="https://rogue.jesed.dev"><img src="https://img.shields.io/badge/site-rogue.jesed.dev-blue" alt="site" /></a>
  <a href="https://jesed.dev"><img src="https://img.shields.io/badge/by-jesed.dev-111" alt="jesed.dev" /></a>
</p>

---

## What it does

Rogue predicts extreme wave events — **rogue ocean waves** (≥ 2.2× significant wave height),
clear-air turbulence, plasma instabilities, and fiber-optic damage — from the
breakthrough mathematics of **finite-time blow-up in supercritical defocusing NLS**
(Merle–Raphael–Rodnianski–Szeftel, 2026 Breakthrough Prize in Mathematics).

Two surfaces: a **browser dashboard** (no install) and a **Rust CLI** (for servers,
automation, streaming instruments). Both run the same symplectic split-step
NLS solver — deterministic, seeded, and reproducible.

## Why

- **$2B+ / year** in shipping damage from rogue waves
- **~100 cruise-ship deaths / year**
- Offshore platform, submarine cable, and wind turbine damage
- Existing tools (WAVEWATCH III, SWAN) are spectral — they don't resolve individual waves
- CFD is accurate but far too slow for real-time
- Statistical methods (Rayleigh) systematically underestimate extremes

Finite-time blow-up detection gives a mathematically rigorous, real-time flag.

---

- [Live site] &nbsp;·&nbsp; [Education dashboard] &nbsp;·&nbsp; [Production forecast]
- [Download binary] &nbsp;·&nbsp; [crates.io] &nbsp;·&nbsp; [Contributing]

---

## Quickstart

### No Rust — download a binary

Pre-built binaries for Linux, macOS, and Windows are attached to every
[GitHub Release]. Download, extract, and run:

```bash
chmod +x rogue
./rogue predict data/sample-buoy.csv
./rogue help
```

### Web (browser, no build)

```bash
cd education && npm install && npm run dev
```

Open [localhost:5173](http://localhost:5173) — the **home page** links to the
interactive education playground and the production forecast form.

### From source (Rust)

```bash
git clone https://github.com/jesedv/rogue.git && cd rogue
./dev setup     # cargo fetch + npm install + wasm-pack build
./dev test      # cargo test --workspace
./dev ui        # browser dashboard at localhost:5173
```

## CLI

```bash
./dev b                           # release build
./dev r predict sea.csv           # production forecast
./dev r predict sea.csv --json    # machine-readable JSON
./dev r akhmediev 1.0 1.0         # Akhmediev breather
./dev r peregrine 1.0             # Peregrine soliton
./dev r ocean 60.0                 # JONSWAP random sea
./dev r blowup 7.0 0.6            # supercritical blow-up
./dev r soliton 4.0               # soliton walk
```

Or install system-wide:

```bash
cargo install --path . --root ~/.local
rogue predict data/sample-buoy.csv --json
```

## Production — live data intake

`rogue predict` reads a CSV/TSV of sea-state observations (buoys, AIS, bridge,
streaming instruments) and outputs a physical, dimensionality-correct forecast
per row.

**Input format** (`t, hs, tp` — optional `gamma`):

```csv
t,hs,tp,gamma
0,1.5,8.0,3.3
3600,2.8,9.5,3.3
```

Column aliases accepted: `time`, `timestamp`, `Hs`, `significant_wave_height`,
`Tp`, `peakperiod`, `gamma`.

```bash
./dev p data/sample-buoy.csv    # human-readable report
./dev p data/sample-buoy.csv --json  # JSON for alarms / routing systems
```

The bridge derives carrier wavenumber `k₀ = ω₀²/g` (deep-water), group-velocity
dispersion `β = ω₀/(8k₀²)` m²/s, and nonlinearity `γ = ω₀k₀²/2` /m·s from
`Tp`. Output is dimensional: crest factor in σ, significant wave height in
metres.

## Connecting instruments

Any buoy, bridge, AIS feed, or streaming source that writes `t, hs, tp` can
feed the solver directly:

```bash
watch -n 60 'rogue predict /var/feed/current.csv --json >> alerts.log'
```

The forecast is deterministic per seed — wind back an AIS log and the alerts
reproduce identically for audit.

## Repository layout

```
crates/
  rogue-nls          # split-step NLS solver (1D/2D), FFT, diagnostics, scenarios
  rogue-blow-up      # virial / energy + H¹-acceleration blow-up detection, ETA
  rogue-rogue        # Akhmediev / KM / Peregrine breathers, Benjamin–Feir, crest detection
  rogue-fluid        # Madelung ↔ compressible Euler bridge
  rogue-plasma       # plasma modulational instability
  rogue-fiber        # fiber-optic soliton power / bandwidth (supercontinuum)
  rogue-turbulence   # clear-air turbulence intermittency
  rogue-wasm         # wasm-bindgen bindings for the browser dashboard
  rogue-production   # real-data intake + physical forecast bridge
education/           # interactive browser dashboard (Vite + TS + WASM)
data/                # sample sea-state CSV records
docs/math.md         # full mathematical basis
scripts/             # regression + bench suite (including production regress)
```

## Regression suite

```bash
./dev regress          # Akhmediev / Peregrine / ocean / blowup physics suite
./dev regress-prod     # production CSV-intake regression
```

## Benchmarks (release, 12 threads)

| Constraint | Target | Measured |
|---|---|---|
| 2D step 1024×1024 | < 100 ms | 59 ms |
| Blow-up detection | < 1 s | streaming fit |
| WASM bundle | ≤ 8 MB | verified at build |

## Math

`docs/math.md` — governing equation `i∂tψ = −(β₂/2)∇²ψ + γ|ψ|^{p−1}ψ`,
invariants, Strang splitting, blow-up detection (virial + H¹ ETA),
Benjamin–Feir law `γ(ν) = ν√(2A² − ν²)`, unified breather formula, Madelung /
Euler bridge.

## License

**Free software.** MIT — do whatever you want. Built for maritime, aviation,
plasma, and optics operators worldwide with zero restrictions.

## Contributing

Contributions are welcome. Open an issue, fork the repo, and send a PR.

- Branch from `main`
- Run `./dev check` before pushing (cargo check + test)
- For physics changes, include a regression run (`./dev regress`)
- For production changes, include `./dev regress-prod`

## Links

- **Live site:** [rogue.jesed.dev](https://rogue.jesed.dev)
- **Author:** [jesed.dev](https://jesed.dev)
- **GitHub:** [github.com/jesedv/rogue](https://github.com/jesedv/rogue)
- **crates.io:** [rogue-nls](https://crates.io/crates/rogue-nls)
- **Docs:** `docs/math.md`

[Live site]: https://rogue.jesed.dev
[Education dashboard]: https://rogue.jesed.dev/education
[Production forecast]: https://rogue.jesed.dev/production
[Download binary]: https://github.com/jesedv/rogue/releases
[GitHub Release]: https://github.com/jesedv/rogue/releases
[crates.io]: https://crates.io/crates/rogue-nls
[Contributing]: #contributing