"""
    Constants

Physical constants and reference state for line-by-line absorption. Base values are
defined once in `constants.toml` and bound here to concrete `Float64` globals for
type-stable, allocation-free access; derived mathematical constants are computed
from them. Profile math converts to the working `FT` at the call site, so the
Float32/Metal paths never promote.
"""
module Constants

using TOML

export LN2, SQRT_LN2, SQRT_LN2_OVER_SQRT_PI, SQRT_2LN2, SQRT_PI,
       C2_RAD, K_BOLTZMANN, C_LIGHT, AMU, P_REF, T_REF, NM_PER_M

const _TOML_PATH = joinpath(@__DIR__, "constants.toml")
include_dependency(_TOML_PATH)              # recompile when the definitions change
const _DEF = TOML.parsefile(_TOML_PATH)

# Float64(...) (not a ::Float64 assertion) so an integer-valued TOML entry like
# `temperature = 296` still loads instead of throwing a TypeError.
"Second radiation constant hc/k_B [cm·K]."
const C2_RAD      = Float64(_DEF["physical"]["c2_radiation"])
"Boltzmann constant [J/K]."
const K_BOLTZMANN = Float64(_DEF["physical"]["k_boltzmann"])
"Speed of light [m/s]."
const C_LIGHT     = Float64(_DEF["physical"]["c_light"])
"Atomic mass unit [kg]."
const AMU         = Float64(_DEF["physical"]["amu"])
"Reference pressure [hPa]."
const P_REF       = Float64(_DEF["reference"]["pressure"])
"Reference temperature [K]."
const T_REF       = Float64(_DEF["reference"]["temperature"])
"Wavelength[nm] = `NM_PER_M` / wavenumber[cm⁻¹]."
const NM_PER_M    = Float64(_DEF["conversion"]["nm_per_m"])

"Natural log of 2."
const LN2                   = log(2)
"√ln2."
const SQRT_LN2              = sqrt(LN2)
"√(2 ln2) — relates Doppler HWHM to the Gaussian standard deviation."
const SQRT_2LN2             = sqrt(2 * LN2)
"√ln2 / √π — the Doppler/Voigt area normalization."
const SQRT_LN2_OVER_SQRT_PI = SQRT_LN2 / sqrt(π)
"√π."
const SQRT_PI               = sqrt(π)

end # module
