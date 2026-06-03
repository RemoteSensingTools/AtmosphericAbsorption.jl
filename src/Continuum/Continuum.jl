"""
    Continuum

Continuum absorption: HITRAN collision-induced absorption (CIA) and the MT_CKD water-
vapor continuum. Both expose cross-section quantities (CIA: binary σ_AB(ν,T)
[cm⁵/molec²]; MT_CKD: σ(ν) [cm²/molec]); combining them with number densities and path
length to get optical depth is left to the radiative-transfer caller.
"""
module Continuum

using Arrow
using ..Constants: C2_RAD

export CIABlock, CIATable, parse_cia_file, build_cia_table, load_cia,
       cia_cross_section, cia_cross_section!,
       MTCKDTable, MTCKDBand, load_mtckd, build_mtckd_band, h2o_continuum, h2o_continuum!

include("cia.jl")
include("mtckd.jl")

end # module
