# Line shapes

The shape of a spectral line — how absorption is distributed in wavenumber around a transition's center ν₀ — is set by the physics of the gas: how fast the molecules move (thermal Doppler broadening) and how often they collide (pressure broadening). AtmosphericAbsorption.jl lets you select the broadening physics with a single `profile=` keyword on `LineByLineModel`, while keeping the same line list, partition function, and grid. This page explains the profile family and how to choose among them.

<iframe title="Line-shape families" src="assets/plots/lineshape_families.html" loading="lazy" style="width:100%;height:520px;border:1px solid var(--vp-c-divider);border-radius:8px;"></iframe>

## The profile family

Each profile is a more complete physical model than the one before it. Pass any of them as `profile=...`:

| Profile | Constructor | Physics captured |
|---------|-------------|------------------|
| Doppler | `Doppler()` | Thermal motion only (zero-pressure limit). Gaussian. |
| Lorentz | `Lorentz()` | Collisional (pressure) broadening only. Lorentzian. |
| Voigt | `Voigt()` | Convolution of Doppler and Lorentz — the standard workhorse. |
| Speed-dependent Voigt | `SpeedDependentVoigt()` | Voigt plus the speed dependence of collisional width/shift. |
| Hartmann-Tran | `HartmannTran()` | Full pCqSDHC: speed dependence + velocity-changing (Dicke) collisions. |

### Doppler — thermal broadening

At vanishing pressure the only broadening mechanism is the thermal velocity distribution of the absorbers. The result is a Gaussian whose half-width γ_d scales with √(T/M). `Doppler()` is the zero-pressure limit of the Voigt profile and is mostly useful for upper-atmosphere / line-core work and teaching.

### Lorentz — pressure broadening

Collisions interrupt the radiating dipole and broaden the line into a Lorentzian with a half-width γ_l that grows with pressure. This is the opposite limit — pure collisional broadening with no thermal contribution. In the troposphere this dominates the wings.

### Voigt — the convolution

Real lines are both Doppler- and pressure-broadened, so the observed shape is the convolution of a Gaussian and a Lorentzian: the **Voigt profile**. This is the standard choice for almost all atmospheric work and is the default:

```julia
model = LineByLineModel(lines, partition; profile=Voigt())
```

The convolution has no closed form; it is evaluated through the **complex probability function** (CPF) — see [CPF strategy](#cpf-strategy-advanced) below.

### Speed-dependent Voigt

The collisional width and shift actually depend on the relative speed of the colliding molecules, which the plain Voigt model ignores. Accounting for this speed dependence narrows and subtly reshapes the core. Select it with `SpeedDependentVoigt()`. The advanced line parameters it needs (γ₂, δ₂) come from `load_hitran_nonvoigt`.

### Hartmann-Tran — the full pCqSDHC model

The Hartmann-Tran profile (the partially-correlated quadratic speed-dependent hard-collision model, **pCqSDHC**) adds velocity-changing (Dicke narrowing) collisions on top of the speed dependence. It is the most complete line shape in routine use and is recommended when matching reference codes at the highest accuracy. Select it with `HartmannTran()`. It consumes the full set of advanced columns (γ₂, δ₂, η, ν_VC).

## Getting the advanced parameters

`Doppler()`, `Lorentz()`, and `Voigt()` run from a standard HITRAN or ExoMol line list. `SpeedDependentVoigt()` and `HartmannTran()` need the extra Hartmann-Tran columns, which come from HITRAN's authenticated non-Voigt API:

```julia
activate_hitran!("your-key")   # or set ENV["HITRAN_API_KEY"]
db = load_hitran_nonvoigt("H2O"; numin=3700, numax=3850, min_strength=0.0, FT=Float64)
```

This fills the database's advanced columns (`γ2_air`, `δ2_air`, `η`, `νVC`), after which `SpeedDependentVoigt()` and `HartmannTran()` have everything they need.

## Line mixing folds in automatically

Voigt, Lorentz, speed-dependent Voigt, and Hartmann-Tran are **collisional** profiles, and they apply **first-order (Rosenkranz) line mixing** automatically using each line's Y_LM coefficient. Line mixing matters where neighboring lines overlap strongly (CO₂ Q-branches, O₂ A-band), redistributing intensity between the lines and reshaping the band envelope. You do not enable it separately — it is built into the collisional profiles. `Doppler()`, being the zero-pressure limit, carries no line mixing.

## Swapping profiles on the same line list

Because the profile is just a keyword, you can hold the line list, partition function, grid, pressure, and temperature fixed and swap only the broadening physics:

```julia
using AtmosphericAbsorption

# Authenticated load so the advanced (HT) columns are populated
activate_hitran!("your-key")
lines     = load_hitran_nonvoigt("H2O"; numin=3700, numax=3850, FT=Float64)
partition = TIPS2017PF()

grid = collect(3700.0:0.01:3850.0)   # cm⁻¹
p, T = 1013.25, 296.0                 # hPa, K

# Same line list, three different physics models
σ_voigt = compute_cross_section(
    LineByLineModel(lines, partition; profile=Voigt()),               grid, p, T)
σ_sdv   = compute_cross_section(
    LineByLineModel(lines, partition; profile=SpeedDependentVoigt()), grid, p, T)
σ_ht    = compute_cross_section(
    LineByLineModel(lines, partition; profile=HartmannTran()),        grid, p, T)
```

Each `σ` is a `Vector` in cm²/molecule on the model's device. The differences between `σ_voigt`, `σ_sdv`, and `σ_ht` are exactly the speed-dependence and Dicke-narrowing effects described above.

## CPF strategy (advanced)

The Voigt-family profiles are all evaluated through a **complex probability function** (the Faddeeva function w(z)). The CPF strategy is selectable but rarely needs changing:

- `HumlicekWeideman32()` — default. A 32-term rational approximation that is fast and GPU-safe.
- `ErfcxCPF()` — a CPU reference implementation using `erfcx`, useful for cross-checking accuracy.

```julia
model = LineByLineModel(lines, partition; profile=Voigt(), cpf=ErfcxCPF())
```

## Low-level shape functions (teaching)

For plots and pedagogy you can call the bare shape functions directly, without building a model. They take detunings Δν = ν − ν₀ and the relevant widths (all in cm⁻¹):

```julia
using AtmosphericAbsorption

Δν = -0.5:0.001:0.5

g_doppler = doppler.(Δν, 0.05)            # Gaussian, half-width γ_d
g_lorentz = lorentz.(Δν, 0.05)            # Lorentzian, half-width γ_l

cpf = HumlicekWeideman32()
y   = 0.6                                  # γ_l / γ_d (the Voigt y parameter)
g_voigt = voigt.(Ref(cpf), Δν, 0.05, y)   # Voigt via the CPF

# Full pCqSDHC kernel (returns complex; the real part is the profile)
z = pcqsdhc.(Ref(cpf), 0.0, 0.05, 0.01, 0.0, 0.0, 0.0, 0.0, Δν)
g_ht = real.(z)
```

`pcqsdhc(cpf, ν0, Γ0, Γ2, Δ0, Δ2, νVC, η, ν)` is the single kernel underneath both `SpeedDependentVoigt()` (set νVC = η = 0) and `HartmannTran()` (full arguments). Setting Γ₂ = Δ₂ = νVC = η = 0 reduces it to a plain Voigt, which is a good way to see how the family nests.

## Validation

The line-shape implementations are validated against **HAPI**, the HITRAN reference Python interface:

- Voigt cross-sections on real CO₂/H₂O match HAPI to ≤ 5×10⁻³.
- Hartmann-Tran and speed-dependent Voigt match HAPI's `pcqsdhc` to ~1×10⁻⁶.
- GPU results equal CPU results to machine precision, and the Float32 pipeline matches Float64 to machine precision.

So you can move up the profile family — Voigt → speed-dependent Voigt → Hartmann-Tran — purely on physical-accuracy grounds, confident that each shape reproduces the HITRAN reference.
