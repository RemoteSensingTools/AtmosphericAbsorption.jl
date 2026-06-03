#=
HITRAN line-list backend: parse a `.par` file into the uniform `LineDatabase`,
fill per-line molar mass from the isotopologue table, sort by wavenumber, and report
TIPS-2017 as the partition function.
=#

"""
    HitranPort(path; edition="HITRAN2016")

A HITRAN `.par` file as a line-list source. Advanced (HT/SDV/line-mixing) columns
are left at zero — HITRAN `.par` is Voigt-parameterized.
"""
struct HitranPort <: AbstractLineListPort
    path::String
    edition::String
end

HitranPort(path::AbstractString; edition::AbstractString = "HITRAN2016") =
    HitranPort(String(path), String(edition))

"""
    load_lines(port::HitranPort; mol=-1, iso=-1, ν_min=0.0, ν_max=Inf,
               min_strength=0.0, FT=Float64) -> LineDatabase{FT}

Parse and filter the port's `.par` file into a `LineDatabase{FT}` sorted ascending in
wavenumber, with per-line molar mass attached from the isotopologue table.
"""
function load_lines(port::HitranPort; mol::Integer = -1, iso::Integer = -1,
                    ν_min::Real = 0.0, ν_max::Real = Inf, min_strength::Real = 0.0,
                    FT::Type{<:AbstractFloat} = Float64)
    c = parse_par(port.path; mol, iso, ν_min, ν_max, min_strength)
    p = sortperm(c.ν)
    mm = FT[molar_mass(c.mol[i], c.iso[i]) for i in p]
    return LineDatabase(; mol = c.mol[p], iso = c.iso[p], ν0 = FT.(c.ν[p]), S = FT.(c.S[p]),
                        E_lower = FT.(c.E[p]), g_upper = FT.(c.g_up[p]),
                        γ_air = FT.(c.γ_air[p]), γ_self = FT.(c.γ_self[p]),
                        n_air = FT.(c.n_air[p]), δ_air = FT.(c.δ[p]), molar_mass = mm,
                        meta = source_metadata(port, mol, iso))
end

"""
    partition_function(::HitranPort, mol, iso) -> TIPS2017PF

HITRAN uses TIPS-2017 partition sums (iso-aware via `Q_ratio(pf, mol, iso, …)`).
"""
partition_function(::HitranPort, mol::Integer, iso::Integer) = TIPS2017PF()

"""
    source_metadata(port::HitranPort, mol, iso) -> SourceMetadata

Reference state (296 K, 1013.25 hPa) and provenance string for the port.
"""
source_metadata(port::HitranPort, mol::Integer, iso::Integer) =
    SourceMetadata("HITRAN $(port.edition)", T_REF, P_REF)
