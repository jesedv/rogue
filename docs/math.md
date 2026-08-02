# Rogue — the mathematics

Finite-time blow-up detection for supercritical *defocusing* NLS
(Merle–Raphael–Rodnianski–Szeftel 2026) applied to rogue-wave prediction.

## 1. Governing equations

The solver integrates the NLS in the physical form

```text
i ∂t ψ = −(β2/2) ∇²ψ + γ |ψ|^(p−1) ψ,        γ = −|γ| focusing, +|γ| defocusing
```

| Name | β2 | γ | power p | Equation |
|---|---|---|---|---|
| `defocusing()` | 2 | +1 | p | `i ∂t ψ + ∇²ψ − |ψ|^(p−1)ψ = 0` (MRRS) |
| `focusing()` | 2 | −1 | p | `i ∂t ψ + ∇²ψ + |ψ|^(p−1)ψ = 0` |
| `fiber(b2,g,…)` | β2 | g | 3 | fiber NLSE |

The MRRS convention makes the defocusing blow-up transparent: for
supercritical power `p > 1 + 4/d` the `L²`-supercritical regime admits
smooth solutions that nonetheless develop a **finite-time singularity**
(the 2026 breakthrough), detected here via virial dynamics and `H¹`
acceleration.

## 2. Invariants (reported by the solver)

For cubic NLS in 1D, with `M = ∫|ψ|² dx`, `E = ∫ (β2/2|ψ_x|² − γ/(p+1)·|ψ|^(p+1)) dx`:

- **Mass** `M` — conserved by construction (unitary split-step).
- **Momentum** `P = ∫ Im(conj ψ · ψ_x) dx` — conserved.
- **Energy** `E` — conserved to `O(dt²)` by the symplectic Strang split.
- **H¹ seminorm** `‖∇ψ‖_L²` — the blow-up monitor (grows without bound).

Parseval scalings used by `Diagnostics` (standard DFT bin order, `dx = lx/nx`):

```text
M = dx Σ |ψ|²
P = (dx/nx) Σ k |ψ̂|²
‖∇ψ‖² = (lx/nx²) Σ k² |ψ̂|²
```

## 3. Time integration — Strang splitting (symplectic)

Each step: linear half-step `dt/2` → nonlinear half-step `dt` → linear
half-step `dt/2`. The linear propagator is exact in Fourier space
(`exp(−i(β2/2)k² dt/2)`), the nonlinear one is a pure phase
(`exp(−iγ|ψ|^(p−1) dt)`). Because both factors are phase rotations, the map is
unitary ⇒ mass exactly conserved, energy conserved to `O(dt²)`. This is the
*geometric integration* requirement for long-time (and pre-blow-up) fidelity
(Hairer–Lubich–Wanner).

Wavenumbers follow DFT bin order: `k_j = 2π·m_j/lx`, `m_j = j` for
`j ≤ n/2`, `m_j = j − n` otherwise (positive band first half, negative second).

## 4. Finite-time blow-up detection (`rogue-blow-up`)

### 4.1 Rigorous (energy) criterion

For focusing NLS with `E < 0` (or `E = 0` with nonzero momentum), blow-up is
guaranteed: the virial identity

```text
d²/dt² ∫ x²|ψ|² dx = 4E + 4(d − (d+2)(p−1)/2)·∫|ψ|² dx · …   (focusing supercritical)
```

forces a focusing collapse. The detector flags `guaranteed_by_energy`
immediately.

### 4.2 Empirical (H¹ acceleration) criterion

Near blow-up the solution concentrates:

```text
‖∇ψ‖(t)² ~ 1/(T* − t)      ⟹      ETA = −a/b  from  y(t) = 1/‖∇ψ‖² ≈ a + b·t
```

A least-squares fit of `y(t) = 1/‖∇ψ‖²` over a sliding window gives the
blow-up time `T*`. `active` when energy is negative **or** the log-growth
rate of `H¹` exceeds threshold with a valid ETA fit.

## 5. Rogue waves (`rogue-rogue`)

### 5.1 Benjamin–Feir / modulational instability

Solver convention `iψ_t + ψ_xx + |ψ|²ψ = 0`, background amplitude `A`,
sideband wavenumber `ν`:

```text
γ(ν) = ν · √(2A² − ν²),   unstable for 0 < ν < √2·A,  max γ = A² at ν = A
```

Verified numerically against the split-step solver across the instability band
(rel. err < 1.5%, `examples/probe_bf.rs`).

### 5.2 Exact breathers

Reference convention `iu_t + u_xx + 2|u|²u = 0` (unit background `e^{2it}`);
mapped to the solver by `ψ(x,t) = √2·u(x,t)` (amplitude rescale, *not* a time
rescale). Unified breather (`p = 2 sin φ`, `Ω = 2 sin 2φ`):

```text
u(x,t) = [ cosh(Ωt − 2iφ) − cos(φ)cos(px) ] / [ cosh(Ωt) − cos(φ)cos(px) ] · e^{2it}
```

| φ | Family | Behavior |
|---|---|---|
| real, `φ ∈ (0, π/2)` | Akhmediev | periodic in `x` (period `π/sin φ`), breathes once; crest `≤ 3×` |
| `φ = iφ₀` | Kuznetsov–Ma | localized in `x`, periodic in `t` |
| `φ → 0` | Peregrine | `u = e^{2it}(1 − 4(1+4it)/(1+4x²+16t²))`, crest exactly `3×` |

`akhmmediev_field`/`peregrine_field` produce solver-valid initial conditions;
`solver_reproduces_akhmediev_breather` replays the breather end-to-end.

### 5.3 Crest criterion

A rogue event is recorded when the local crest exceeds `RogueDetector`'s
threshold (default `2.2 ×` the significant amplitude `σ = rms(|ψ|)`), matching
the oceanographic `H_max ≥ 2.2·H_s` definition.

## 6. Fluid bridge (`rogue-fluid`)

Madelung transform `ψ = √ρ e^{iφ}`, `u = ∂xφ`. The defocusing cubic NLS is
isothermal compressible Euler plus quantum pressure:

```text
ρ_t + (ρu)_x = 0
u_t + uu_x + ρ_x = (1/2) ∂x( (√ρ)_xx / √ρ )
```

When the quantum-pressure term is negligible this is exactly isothermal Euler
(sound speed 1); `madelung()` extracts `(ρ, u, q)` from any field and
`euler.rs` provides a Rusanov–Heun Euler solver cross-checked against NLS.

## 7. Physics platforms

- **`rogue-fiber`** — NLSE soliton power `P = (|β2|/T₀²)/γ`, `N`-soliton energy,
  20 dB spectral bandwidth (supercontinuum / fiber damage).
- **`rogue-plasma`** — modulational growth on a plasma background:
  `γ(ν) = ν√(β2·g·A² − (β2/2)²ν²)`, max `g·A²` at `ν = √(2g/β2)·A`.
- **`rogue-turbulence`** — clear-air turbulence: excess kurtosis of `|ψ|`
  spikes ⇒ intermittent extreme events (severity ladder Smooth→Severe).

## 8. Performance

Hard constraints and current release-build numbers (12 threads):

| Constraint | Target | Measured |
|---|---|---|
| 2D step, 1024×1024 | < 100 ms | 59.3 ms |
| Blow-up ETA detect | < 1 s | trivial (streaming fit) |
| WASM bundle | ≤ 8 MB | verified at build |

## References

- Merle, Raphael, Rodnianski, Szeftel — *Defocusing Wave Singularities* (2026).
- Akhmediev, Ankiewicz, Taki — *Waves that Appear from Nowhere…* Phys. Lett. A 2009.
- Onorato et al. — *Rogue Waves and Their Generating Mechanisms*, Phys. Rep. 2013.
- Benjamin, Feir — *The Disintegration of Wave Trains on Deep Water*, JFM 1967.
- Hairer, Lubich, Wanner — *Geometric Numerical Integration*, Springer 2006.
- Tao — *Nonlinear Dispersive Equations*, CBMS 2006.
