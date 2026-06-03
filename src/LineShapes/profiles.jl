#=
Area-normalized line profiles (∫ ϕ dν = 1) as functions of detuning Δν = ν − ν₀.
The bare `doppler`/`lorentz`/`voigt` are the implementations; `evaluate` gives a
uniform, compile-time-dispatched signature for the cross-section kernel to call per
line. All functions are FT-generic and GPU-safe (constants converted to FT inline).
=#

"""Singleton tag for a line-shape model, dispatched in `evaluate`."""
abstract type AbstractLineProfile end

"""Gaussian thermal (Doppler) profile."""
struct Doppler <: AbstractLineProfile end

"""Lorentzian pressure-broadened profile."""
struct Lorentz <: AbstractLineProfile end

"""Voigt profile — Doppler⊗Lorentz convolution via the complex probability function."""
struct Voigt <: AbstractLineProfile end

"""Speed-dependent Voigt profile (quadratic speed dependence of width/shift)."""
struct SpeedDependentVoigt <: AbstractLineProfile end

"""Hartmann-Tran (HTP) profile — speed dependence + velocity-changing collisions."""
struct HartmannTran <: AbstractLineProfile end

"""Gaussian profile of HWHM-related width `γd` (the Doppler HWHM)."""
@inline doppler(Δν::FT, γd::FT) where {FT} =
    FT(SQRT_LN2_OVER_SQRT_PI) / γd * exp(-FT(LN2) * (Δν / γd)^2)

"""Lorentzian profile of HWHM `γl`."""
@inline lorentz(Δν::FT, γl::FT) where {FT} = γl / (FT(π) * (γl^2 + Δν^2))

"""Voigt profile from Doppler HWHM `γd` and width ratio `y = √ln2·γl/γd`."""
@inline voigt(cpf::AbstractCPF, Δν::FT, γd::FT, y::FT) where {FT} =
    FT(SQRT_LN2_OVER_SQRT_PI) / γd * real(w(cpf, FT(SQRT_LN2) / γd * Δν + im * y))

# Uniform per-line entry point: `νI` is the grid wavenumber, `ν0` the line center, and
# `p` a NamedTuple of prepared parameters (γd, Γ0, Γ2, Δ0, Δ2, νVC, η, Y). Each profile reads
# only what it needs; the cross-section kernel calls the same signature for every profile,
# specialized at compile time on the profile type. Δ0 is the pressure shift; `Y` is the
# first-order (Rosenkranz) line-mixing coefficient — it adds the asymmetric `Y·Im` term to
# the collisional profiles and is zero for lists without mixing. Doppler is the zero-
# pressure limit, where mixing vanishes, so it ignores `Y`.
@inline evaluate(::Doppler, ::AbstractCPF, νI::FT, ν0::FT, p) where {FT} =
    doppler(νI - (ν0 + p.Δ0), p.γd)

@inline function evaluate(::Lorentz, ::AbstractCPF, νI::FT, ν0::FT, p) where {FT}
    Δν = νI - (ν0 + p.Δ0)
    return (p.Γ0 + p.Y * Δν) / (FT(π) * (p.Γ0^2 + Δν^2))      # Im part of the complex Lorentzian
end

@inline function evaluate(::Voigt, cpf::AbstractCPF, νI::FT, ν0::FT, p) where {FT}
    Δν = νI - (ν0 + p.Δ0)
    W  = w(cpf, FT(SQRT_LN2) / p.γd * Δν + im * (FT(SQRT_LN2) * p.Γ0 / p.γd))
    return FT(SQRT_LN2_OVER_SQRT_PI) / p.γd * (real(W) + p.Y * imag(W))
end
