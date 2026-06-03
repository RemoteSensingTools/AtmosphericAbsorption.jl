#=
Authenticated HITRANonline (HAPI2) fetch of NON-Voigt line-shape parameters
(Hartmann-Tran / speed-dependent Voigt) and parsing into the LineDatabase advanced
columns that the `pcqsdhc` HT/SDV profiles consume. Two-step API:
  1. GET /api/v2/{key}/transitions?...&request_params=par_line,<HT params>  → JSON header
     naming a data file.
  2. GET /results/{file}  → 160-char par_line + comma-separated requested params per line
     ('#' = parameter absent for that line).
Requires an API key (`activate_hitran!` or `HITRAN_API_KEY`). The key-bearing URL is never
written to disk; provenance records the query, data file, hash and date only.
=#

const _HITRAN_API     = "https://hitran.org/api/v2"
const _HITRAN_RESULTS = "https://hitran.org/results"
# Air-broadened HT parameters at 296 K, in request order → appended-column order.
const _NONVOIGT_PARAMS = ["gamma_HT_0_air_296", "gamma_HT_2_air_296", "delta_HT_0_air_296",
                          "delta_HT_2_air_296", "eta_HT_air", "nu_HT_air"]

# Parse an appended parameter token ('#' or empty ⇒ missing, returned as NaN).
@inline function _nv(tokens, k)
    k ≤ length(tokens) || return NaN
    t = strip(tokens[k])
    (isempty(t) || t == "#") ? NaN : parse(Float64, t)
end

"""
    fetch_hitran_nonvoigt(molecule; numin, numax, edition="HITRAN-HT", force=false) -> path

Download HITRAN line data **with non-Voigt (HT/SDV) parameters** for `molecule` over
`[numin, numax]` via the authenticated API, caching the result with provenance. Needs an
API key. Returns the local data-file path.
"""
function fetch_hitran_nonvoigt(molecule::AbstractString; numin::Real = 0, numax::Real = 150000,
                               edition::AbstractString = "HITRAN-HT", force::Bool = false)
    dir = joinpath(_hitran_dir(), edition); mkpath(dir)
    data, meta = joinpath(dir, "$molecule.data"), joinpath(dir, "$molecule.meta.toml")
    # Reuse a covering cache without needing an API key (offline/keyless re-runs).
    if !force && isfile(data) && _hitran_cache_covers(meta, numin, numax)
        return data
    end

    key = hitran_api_key()
    ids = _hitran_iso_ids(molecule_number(molecule))
    isempty(ids) && error("no HITRAN isotopologues found for \"$molecule\"")
    req = join(["par_line"; _NONVOIGT_PARAMS], ",")
    # Step 1: header (key in URL — kept in memory, never stored).
    header_url = "$(_HITRAN_API)/$(key)/transitions?iso_ids_list=$(join(ids, ','))" *
                 "&numin=$(numin)&numax=$(numax)&head=false&fixwidth=0&request_params=$(req)"
    buf = IOBuffer()
    Downloads.request(header_url; output = buf, throw = false)
    m = match(r"\"data\":\s*\"([^\"]+)\"", String(take!(buf)))
    m === nothing &&
        error("HITRAN non-Voigt: no data file in the API response for $molecule " *
              "(check the API key and that the molecule has HT parameters in this window).")

    # Step 2: the (public) results file.
    open(data, "w") do out
        Downloads.request("$(_HITRAN_RESULTS)/$(m.captures[1])"; output = out, throw = false)
    end
    (isfile(data) && filesize(data) > 0 && _looks_like_par(data)) ||
        (rm(data; force = true); error("HITRAN non-Voigt download empty/invalid for $molecule."))

    open(meta, "w") do io          # provenance WITHOUT the key/authenticated URL
        println(io, "molecule = \"", molecule, "\"")
        println(io, "edition = \"", edition, "\"")
        println(io, "numin = ", numin)
        println(io, "numax = ", numax)
        println(io, "request_params = \"", req, "\"")
        println(io, "data_file = \"", m.captures[1], "\"")
        println(io, "sha256 = \"", bytes2hex(SHA.sha256(read(data))), "\"")
        println(io, "fetched = \"", Dates.today(), "\"")
    end
    return data
end

"""
    load_hitran_nonvoigt(molecule; numin=0, numax=Inf, min_strength=0.0, FT=Float64,
                         edition="HITRAN-HT", force=false) -> LineDatabase{FT}

Fetch and parse HITRAN data with non-Voigt parameters into a `LineDatabase{FT}` whose
advanced columns (γ2_air, δ2_air, η, νVC) are populated where HITRAN provides them; the
speed-averaged width/shift (γ_air, δ_air) use the HT values when present, else the Voigt
ones. Pair with `profile = SpeedDependentVoigt()` or `HartmannTran()`.
"""
function load_hitran_nonvoigt(molecule::AbstractString; numin::Real = 0, numax::Real = Inf,
                              min_strength::Real = 0.0, FT::Type{<:AbstractFloat} = Float64,
                              edition::AbstractString = "HITRAN-HT", force::Bool = false)
    path = fetch_hitran_nonvoigt(molecule; numin, numax = isfinite(numax) ? numax : 150000,
                                 edition, force)
    mol = Int32[]; iso = Int32[]; ν0 = Float64[]; S = Float64[]; E = Float64[]; gup = Float64[]
    γair = Float64[]; γself = Float64[]; nair = Float64[]; δair = Float64[]
    γ2 = Float64[]; δ2 = Float64[]; η = Float64[]; νVC = Float64[]; mm = Float64[]
    for line in eachline(path)
        length(line) < 160 && continue
        ν = _f64(line, _PAR.ν)
        numin ≤ ν ≤ numax || continue
        s = _f64(line, _PAR.S)
        s ≥ min_strength || continue
        m = parse(Int32, @view line[_PAR.mol]); i = Int32(_isoid(line[3]))
        t = length(line) > 160 ? split(line[161:end], ',') : SubString{String}[]
        g0, g2, d0, d2, et, nv = _nv(t, 2), _nv(t, 3), _nv(t, 4), _nv(t, 5), _nv(t, 6), _nv(t, 7)
        pγ, pδ = _f64(line, _PAR.γ_air), _f64(line, _PAR.δ)
        push!(mol, m); push!(iso, i); push!(ν0, ν); push!(S, s)
        push!(E, _f64(line, _PAR.E)); push!(gup, _f64(line, _PAR.g_up))
        push!(γair, isnan(g0) ? pγ : g0); push!(γself, _f64(line, _PAR.γ_self))
        push!(nair, _f64(line, _PAR.n_air)); push!(δair, isnan(d0) ? pδ : d0)
        push!(γ2, isnan(g2) ? 0.0 : g2); push!(δ2, isnan(d2) ? 0.0 : d2)
        push!(η, isnan(et) ? 0.0 : et); push!(νVC, isnan(nv) ? 0.0 : nv)
        push!(mm, molar_mass(m, i))
    end
    p = sortperm(ν0)
    return LineDatabase(; mol = mol[p], iso = iso[p], ν0 = FT.(ν0[p]), S = FT.(S[p]),
                        E_lower = FT.(E[p]), g_upper = FT.(gup[p]), γ_air = FT.(γair[p]),
                        γ_self = FT.(γself[p]), n_air = FT.(nair[p]), δ_air = FT.(δair[p]),
                        γ2_air = FT.(γ2[p]), δ2_air = FT.(δ2[p]), η = FT.(η[p]), νVC = FT.(νVC[p]),
                        molar_mass = FT.(mm[p]),
                        meta = SourceMetadata("HITRAN-HT $molecule", T_REF, P_REF))
end
