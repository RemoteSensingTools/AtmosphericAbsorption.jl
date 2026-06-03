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

"""Gaussian profile of HWHM-related width `γd` (the Doppler HWHM)."""
@inline doppler(Δν::FT, γd::FT) where {FT} =
    FT(SQRT_LN2_OVER_SQRT_PI) / γd * exp(-FT(LN2) * (Δν / γd)^2)

"""Lorentzian profile of HWHM `γl`."""
@inline lorentz(Δν::FT, γl::FT) where {FT} = γl / (FT(π) * (γl^2 + Δν^2))

"""Voigt profile from Doppler HWHM `γd` and width ratio `y = √ln2·γl/γd`."""
@inline voigt(cpf::AbstractCPF, Δν::FT, γd::FT, y::FT) where {FT} =
    FT(SQRT_LN2_OVER_SQRT_PI) / γd * real(w(cpf, FT(SQRT_LN2) / γd * Δν + im * y))

# Uniform per-line entry point. Unused widths are ignored at compile time, so the
# cross-section kernel can pass the same argument list for every profile.
@inline evaluate(::Doppler, ::AbstractCPF, Δν::FT, γd::FT, γl::FT, y::FT) where {FT} = doppler(Δν, γd)
@inline evaluate(::Lorentz, ::AbstractCPF, Δν::FT, γd::FT, γl::FT, y::FT) where {FT} = lorentz(Δν, γl)
@inline evaluate(::Voigt, cpf::AbstractCPF, Δν::FT, γd::FT, γl::FT, y::FT) where {FT} = voigt(cpf, Δν, γd, y)
