#=
HITRAN line-list backend: parse a `.par` file into the uniform `LineDatabase`,
fill per-line molar mass from the isotopologue table, sort by wavenumber, and report
TIPS-2017 as the partition function.
=#

"""
    HitranPort(path; edition="HITRAN2016")    # a local `.par` file
    HitranPort(; edition="HITRAN2020")        # the HITRAN database, fetched on demand

A HITRAN line-list source. Two flavours, both consumed the same way by `load_lines`:

  * `HitranPort("co2.par")` wraps a **local `.par` file** — `load_lines` parses and filters it.
  * `HitranPort(; edition)` is a **remote handle** carrying only the edition; the molecule and
    band live on the `load_lines` call (which downloads them, cached). Define the handle once and
    reuse it to pull as many molecules / bands as you like:

    ```julia
    port = HitranPort(; edition="HITRAN2020")
    co2  = load_lines(port; mol=:CO2, ν_min=6300, ν_max=6400)
    h2o  = load_lines(port; mol=:H2O, ν_min=7000, ν_max=7100)
    ```

Advanced (HT/SDV/line-mixing) columns are left at zero — HITRAN `.par` is Voigt-parameterized;
see [`load_hitran_nonvoigt`](@ref) for those. The empty path is the remote sentinel.
"""
struct HitranPort <: AbstractLineListPort
    path::String
    edition::String
end

HitranPort(path::AbstractString; edition::AbstractString = "HITRAN2016") =
    HitranPort(String(path), String(edition))
HitranPort(; edition::AbstractString = "HITRAN2020") = HitranPort("", String(edition))

"""
    load_lines(port::HitranPort; mol=:ALL, iso=:ALL, ν_min=0.0, ν_max=Inf,
               min_strength=0.0, force=false, FT=Float64) -> LineDatabase{FT}

Produce a `LineDatabase{FT}` sorted ascending in wavenumber, with per-line molar mass and the
TIPS-2021 partition function attached. For a **local-file** port the `.par` is parsed and filtered;
for a **remote** port (`HitranPort(; edition)`) the `mol`/band drive a download (cached) of that
molecule over `[ν_min, ν_max]`, which is then parsed.

`mol`/`iso` accept the generic notation — a symbol (`:CO2`, `:main`), a string, an integer id, or
`:ALL` for every molecule/isotopologue (the multi-molecule "see what's in this band" case; only
valid for a local file — a remote download needs a concrete molecule). `force=true` re-downloads
even if cached.
"""
function load_lines(port::HitranPort; mol = -1, iso = -1,
                    ν_min::Real = 0.0, ν_max::Real = Inf, min_strength::Real = 0.0,
                    force::Bool = false, FT::Type{<:AbstractFloat} = Float64)
    m = resolve_molecule(mol); i = resolve_isotopologue(m, iso)
    path = if isempty(port.path)
        m == -1 && throw(ArgumentError(
            "a remote HitranPort needs a molecule — e.g. load_lines(port; mol=:CO2, ν_min=…, ν_max=…)"))
        fetch_hitran(m; numin = ν_min, numax = isfinite(ν_max) ? ν_max : 150000.0,
                     edition = port.edition, force)
    else
        port.path
    end
    c = parse_par(path; mol = m, iso = i, ν_min, ν_max, min_strength)
    p = sortperm(c.ν)
    mm = FT[molar_mass(c.mol[k], c.iso[k]) for k in p]
    return LineDatabase(; mol = c.mol[p], iso = c.iso[p], ν0 = FT.(c.ν[p]), S = FT.(c.S[p]),
                        E_lower = FT.(c.E[p]), g_upper = FT.(c.g_up[p]),
                        γ_air = FT.(c.γ_air[p]), γ_self = FT.(c.γ_self[p]),
                        n_air = FT.(c.n_air[p]), δ_air = FT.(c.δ[p]), molar_mass = mm,
                        meta = source_metadata(port, m, i), partition = partition_function(port, m, i))
end

"""
    partition_function(::HitranPort, mol, iso) -> TIPS2021PF

HITRAN partition sums from the latest edition, TIPS-2021 (iso-aware via
`Q_ratio(pf, mol, iso, …)`). Use `TIPS2017PF()` explicitly to reproduce the previous
edition (e.g. for HAPI cross-validation).
"""
partition_function(::HitranPort, mol::Integer, iso::Integer) = TIPS2021PF()

"""
    source_metadata(port::HitranPort, mol, iso) -> SourceMetadata

Reference state (296 K, 1013.25 hPa) and provenance string for the port.
"""
source_metadata(port::HitranPort, mol::Integer, iso::Integer) =
    SourceMetadata("HITRAN $(port.edition)", T_REF, P_REF)
