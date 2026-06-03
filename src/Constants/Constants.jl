"""
    Constants

Physical constants and reference state for line-by-line absorption.
Values are CODATA-2018. Constants are stored in `Float64`; profile math converts
to the working `FT` at the call site so Float32/Metal paths never promote.
"""
module Constants

export LN2, SQRT_LN2, SQRT_LN2_OVER_SQRT_PI, SQRT_2LN2,
       C2_RAD, K_BOLTZMANN, C_LIGHT, AMU, P_REF, T_REF, NM_PER_M

const LN2                  = log(2)                 # natural log of 2
const SQRT_LN2             = sqrt(LN2)
const SQRT_2LN2            = sqrt(2 * LN2)
const SQRT_LN2_OVER_SQRT_PI = SQRT_LN2 / sqrt(π)    # Voigt/Doppler normalization

const C2_RAD      = 1.4387768775039337   # second radiation constant hc/k_B [cm·K]
const K_BOLTZMANN = 1.380649e-23         # Boltzmann constant [J/K]
const C_LIGHT     = 2.99792458e8         # speed of light [m/s]
const AMU         = 1.66053906660e-27    # atomic mass unit [kg]

const P_REF    = 1013.25   # reference pressure [hPa]
const T_REF    = 296.0     # reference temperature [K]
const NM_PER_M = 1.0e7     # wavelength[nm] = NM_PER_M / wavenumber[cm⁻¹]

end # module
