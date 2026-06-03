#=
ExoMol line-list backend. ExoMol stores energy states (.states), transitions with
Einstein-A coefficients (.trans), a partition function (.pf), and broadening (.broad).
We compute HITRAN-style line intensities from A + states + partition (verified against
HITRAN CO to <0.05%), assign air broadening by branch, and emit the uniform
`LineDatabase`. ExoMol intensities are per-isotopologue, so we apply the terrestrial
isotopic abundance to match HITRAN conventions. ExoMol carries no pressure shift, so
`δ_air` is left at zero.
=#

const _BAR_PER_ATM = 1.01325   # ExoMol γ is cm⁻¹/bar; LineDatabase γ_air is cm⁻¹/atm

"""
    ExoMolPort(molecule, iso_slug, linelist; hitran_mol, hitran_iso)

An ExoMol line list as a data source, e.g.
`ExoMolPort("CO", "12C-16O", "Li2015"; hitran_mol=5, hitran_iso=1)`. The HITRAN
`(mol, iso)` ids supply the terrestrial abundance and label the output lines.
"""
struct ExoMolPort <: AbstractLineListPort
    molecule::String
    iso_slug::String
    linelist::String
    hitran_mol::Int
    hitran_iso::Int
end

ExoMolPort(molecule, iso_slug, linelist; hitran_mol, hitran_iso) =
    ExoMolPort(String(molecule), String(iso_slug), String(linelist), hitran_mol, hitran_iso)

_exomol_meta(port, version) =
    SourceMetadata("ExoMol $(port.molecule) $(port.linelist) (v$version)", T_REF, P_REF)

# Branch key for a transition under the file's diet code (J may be half-integer).
@inline function _broad_key(diet, J_lo, J_up)
    diet == "m0" || return J_lo                       # a0 (and fallback): key on lower J
    ΔJ = J_up - J_lo
    return ΔJ == 1 ? J_lo + 1 : ΔJ == -1 ? -J_lo : J_lo
end

"""
    exomol_intensity(g_up, A, ν, E_lo, Q, abund) -> Float64

HITRAN-style line intensity [cm/molecule] at 296 K from the ExoMol Einstein-A
coefficient: `abund·g_up·A/(8πc·ν²)·exp(-c₂E_lo/Tref)·(1-exp(-c₂ν/Tref))/Q`, with `c`
in cm/s. `abund` is the terrestrial isotopic abundance (ExoMol intensities are
per-isotopologue; HITRAN's are abundance-weighted). Verified vs HITRAN CO to <0.05%.
"""
@inline function exomol_intensity(g_up, A, ν, E_lo, Q, abund)
    c_cm, c₂, Tref = 100 * C_LIGHT, C2_RAD, T_REF
    return abund * g_up * A / (8π * c_cm * ν^2) *
           exp(-c₂ * E_lo / Tref) * (-expm1(-c₂ * ν / Tref)) / Q
end

# Accumulate one `.trans` stream into the per-line collectors, filtering by window and
# strength. Factored out so it can be tested offline on a fixture without downloading.
function _stream_trans!(cols, io, st, diet, broad, def, abund, Q296, ν_min, ν_max, min_strength)
    N = length(st.E)
    for ln in eachline(io)
        t = split(ln)
        length(t) < 3 && continue
        up = parse(Int, t[1]); lo = parse(Int, t[2]); A = parse(Float64, t[3])
        ν  = length(t) ≥ 4 ? parse(Float64, t[4]) : st.E[up] - st.E[lo]
        (ν_min ≤ ν ≤ ν_max && ν > 0) || continue
        (1 ≤ up ≤ N && 1 ≤ lo ≤ N)   || continue        # skip refs to absent states
        E_lo, g_up = st.E[lo], st.g[up]
        S = exomol_intensity(g_up, A, ν, E_lo, Q296, abund)
        S ≥ min_strength || continue
        key = _broad_key(diet, st.J[lo], st.J[up])
        γ_bar, nexp = isinteger(key) ? get(broad, Int(key), (def.default_γ, def.default_n)) :
                                       (def.default_γ, def.default_n)
        push!(cols.ν0, ν); push!(cols.S, S); push!(cols.E, E_lo); push!(cols.g, g_up)
        push!(cols.γ, γ_bar * _BAR_PER_ATM); push!(cols.n, nexp)
    end
end

"""
    load_lines(port::ExoMolPort; mol=-1, iso=-1, ν_min=0.0, ν_max=Inf,
               min_strength=0.0, FT=Float64) -> LineDatabase{FT}

Fetch the ExoMol files (only the `.trans` chunks overlapping `[ν_min, ν_max]`), stream
the transitions, compute the abundance-weighted line intensity at 296 K, assign air
broadening, and return a `LineDatabase{FT}` sorted in wavenumber. `mol`/`iso` are
accepted for interface symmetry but ExoMol files are single-isotopologue.
"""
function load_lines(port::ExoMolPort; mol::Integer = -1, iso::Integer = -1,
                    ν_min::Real = 0.0, ν_max::Real = Inf, min_strength::Real = 0.0,
                    FT::Type{<:AbstractFloat} = Float64)
    m   = fetch_exomol_meta(port)
    def = parse_def(m.def)
    _write_provenance(m.dir, port, m.def, def.version)
    st  = read_states(m.states, def.nstates)
    pf  = read_pf(m.pf)
    diet, broad = read_broad(m.broad)
    abund = abundance(port.hitran_mol, port.hitran_iso)
    Q296  = CubicSpline(pf.Q, pf.T)(T_REF)

    cols = (ν0 = Float64[], S = Float64[], E = Float64[], g = Float64[], γ = Float64[], n = Float64[])
    foreach(v -> sizehint!(v, 100_000), cols)
    for tp in fetch_trans(port, def, Float64(ν_min), Float64(ν_max))
        _open_text(tp) do io
            _stream_trans!(cols, io, st, diet, broad, def, abund, Q296, ν_min, ν_max, min_strength)
        end
    end

    p, n = sortperm(cols.ν0), length(cols.ν0)
    z = zeros(FT, n)
    return LineDatabase(; mol = fill(Int32(port.hitran_mol), n), iso = fill(Int32(port.hitran_iso), n),
                        ν0 = FT.(cols.ν0[p]), S = FT.(cols.S[p]), E_lower = FT.(cols.E[p]),
                        g_upper = FT.(cols.g[p]), γ_air = FT.(cols.γ[p]), γ_self = copy(z),
                        n_air = FT.(cols.n[p]), δ_air = copy(z), molar_mass = fill(FT(def.mass), n),
                        meta = _exomol_meta(port, def.version))
end

"""
    partition_function(port::ExoMolPort, mol, iso) -> TabulatedPF

The ExoMol `.pf` partition function (Q(T) spline), used for the temperature correction
instead of TIPS. Fetches only the `.pf` file.
"""
function partition_function(port::ExoMolPort, mol::Integer, iso::Integer)
    pf = read_pf(fetch_pf(port))
    return TabulatedPF(pf.T, pf.Q)
end

"""
    source_metadata(port::ExoMolPort, mol, iso) -> SourceMetadata

Reference state (296 K, 1013.25 hPa) and provenance (incl. dataset version) for the
ExoMol line list.
"""
source_metadata(port::ExoMolPort, mol::Integer, iso::Integer) =
    _exomol_meta(port, parse_def(fetch_exomol_meta(port).def).version)
