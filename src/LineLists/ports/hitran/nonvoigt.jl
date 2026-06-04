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
# Extended HITRAN parameters at 296 K, in request order → appended-column order. After
# the Hartmann-Tran terms come the first-order (Rosenkranz) line-mixing coefficients Y for
# the HT and SDV models (HITRAN ships no Y temperature exponent, so n_Y_LM stays 0,
# matching HAPI's YRosen default), then the self-broadening T-exponent and pressure shift
# (absent from the 160-char .par format).
const _NONVOIGT_PARAMS = ["gamma_HT_0_air_296", "gamma_HT_2_air_296", "delta_HT_0_air_296",
                          "delta_HT_2_air_296", "eta_HT_air", "nu_HT_air",
                          "Y_HT_air_296", "Y_SDV_air_296", "n_self", "delta_self"]

# Parse an appended parameter token ('#' or empty ⇒ missing, returned as NaN).
@inline function _nv(tokens, k)
    k ≤ length(tokens) || return NaN
    t = strip(tokens[k])
    (isempty(t) || t == "#") ? NaN : parse(Float64, t)
end

# Cache is reusable only if it covers the window AND was fetched with the same parameter
# set — otherwise a stale file (e.g. one without the line-mixing columns) is silently reused.
function _nonvoigt_cache_ok(meta_path, numin, numax, req)
    _hitran_cache_covers(meta_path, numin, numax) || return false
    meta = try TOML.parsefile(meta_path) catch; return false end
    get(meta, "request_params", nothing) == req
end

"""
    fetch_hitran_nonvoigt(molecule; numin, numax, edition="HITRAN-HT", force=false) -> path

Download HITRAN line data **with non-Voigt (HT/SDV) parameters** for `molecule` over
`[numin, numax]` via the authenticated API, caching the result with provenance. Needs an
API key. Returns the local data-file path.
"""
function fetch_hitran_nonvoigt(molecule::Union{Integer,Symbol,AbstractString}; numin::Real = 0,
                               numax::Real = 150000, edition::AbstractString = "HITRAN-HT",
                               force::Bool = false)
    id = resolve_molecule(molecule)
    id < 1 && error("fetch_hitran_nonvoigt needs a specific molecule (e.g. :H2O), not :ALL")
    name = String(molecule_symbol(id))
    dir = joinpath(_hitran_dir(), edition); mkpath(dir)
    data, meta = joinpath(dir, "$name.data"), joinpath(dir, "$name.meta.toml")
    req = join(["par_line"; _NONVOIGT_PARAMS], ",")
    # Reuse a covering cache without needing an API key (offline/keyless re-runs),
    # but only if it was fetched with the current parameter set.
    if !force && isfile(data) && _nonvoigt_cache_ok(meta, numin, numax, req)
        return data
    end

    key = hitran_api_key()
    ids = _hitran_iso_ids(id)
    isempty(ids) && error("no HITRAN isotopologues found for $name")
    # Step 1: header (key in URL — kept in memory, never stored).
    header_url = "$(_HITRAN_API)/$(key)/transitions?iso_ids_list=$(join(ids, ','))" *
                 "&numin=$(numin)&numax=$(numax)&head=false&fixwidth=0&request_params=$(req)"
    buf = IOBuffer()
    Downloads.request(header_url; output = buf, throw = false)
    m = match(r"\"data\":\s*\"([^\"]+)\"", String(take!(buf)))
    m === nothing &&
        error("HITRAN non-Voigt: no data file in the API response for $name " *
              "(check the API key and that the molecule has HT parameters in this window).")

    # Step 2: the (public) results file.
    open(data, "w") do out
        Downloads.request("$(_HITRAN_RESULTS)/$(m.captures[1])"; output = out, throw = false)
    end
    (isfile(data) && filesize(data) > 0 && _looks_like_par(data)) ||
        (rm(data; force = true); error("HITRAN non-Voigt download empty/invalid for $name."))

    open(meta, "w") do io          # provenance WITHOUT the key/authenticated URL
        println(io, "molecule = \"", name, "\"")
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
ones. The first-order line-mixing coefficient `Y_LM` is filled from HITRAN's `Y_HT_air_296`
(falling back to `Y_SDV_air_296`) where present, so line mixing is active for these lines.
The self-broadening temperature exponent (`n_self`) and pressure shift (`δ_self`) — absent
from the basic `.par` format — are also pulled, defaulting to `n_air` and 0 respectively
when HITRAN does not provide them. Pair with `profile = SpeedDependentVoigt()` or
`HartmannTran()`.
"""
function load_hitran_nonvoigt(molecule::Union{Integer,Symbol,AbstractString}; numin::Real = 0,
                              numax::Real = Inf, min_strength::Real = 0.0,
                              FT::Type{<:AbstractFloat} = Float64,
                              edition::AbstractString = "HITRAN-HT", force::Bool = false)
    name = String(molecule_symbol(resolve_molecule(molecule)))
    path = fetch_hitran_nonvoigt(molecule; numin, numax = isfinite(numax) ? numax : 150000,
                                 edition, force)
    return _parse_nonvoigt_data(path; numin, numax, min_strength, FT,
                                source = "HITRAN-HT $name")
end

# Parse a HITRANonline results file (160-char par_line + the comma-separated
# `_NONVOIGT_PARAMS` per line) into a LineDatabase. Split out from the fetch so it can be
# exercised offline on a bundled golden file.
function _parse_nonvoigt_data(path::AbstractString; numin::Real = 0, numax::Real = Inf,
                              min_strength::Real = 0.0, FT::Type{<:AbstractFloat} = Float64,
                              source::AbstractString = "HITRAN-HT")
    mol = Int32[]; iso = Int32[]; ν0 = Float64[]; S = Float64[]; E = Float64[]; gup = Float64[]
    γair = Float64[]; γself = Float64[]; nair = Float64[]; δair = Float64[]
    nself = Float64[]; δself = Float64[]
    γ2 = Float64[]; δ2 = Float64[]; η = Float64[]; νVC = Float64[]; Y = Float64[]; mm = Float64[]
    for line in eachline(path)
        length(line) < 160 && continue
        ν = _f64(line, _PAR.ν)
        numin ≤ ν ≤ numax || continue
        s = _f64(line, _PAR.S)
        s ≥ min_strength || continue
        m = parse(Int32, @view line[_PAR.mol]); i = Int32(_isoid(line[3]))
        t = length(line) > 160 ? split(line[161:end], ',') : SubString{String}[]
        g0, g2, d0, d2, et, nv = _nv(t, 2), _nv(t, 3), _nv(t, 4), _nv(t, 5), _nv(t, 6), _nv(t, 7)
        yht, ysdv = _nv(t, 8), _nv(t, 9)        # line mixing: prefer HT, fall back to SDV
        y = !isnan(yht) ? yht : (isnan(ysdv) ? 0.0 : ysdv)
        ns, ds = _nv(t, 10), _nv(t, 11)         # self-broadening T-exponent and shift
        nairv = _f64(line, _PAR.n_air); pγ, pδ = _f64(line, _PAR.γ_air), _f64(line, _PAR.δ)
        push!(mol, m); push!(iso, i); push!(ν0, ν); push!(S, s)
        push!(E, _f64(line, _PAR.E)); push!(gup, _f64(line, _PAR.g_up))
        push!(γair, isnan(g0) ? pγ : g0); push!(γself, _f64(line, _PAR.γ_self))
        push!(nair, nairv); push!(δair, isnan(d0) ? pδ : d0)
        push!(nself, isnan(ns) ? nairv : ns); push!(δself, isnan(ds) ? 0.0 : ds)
        push!(γ2, isnan(g2) ? 0.0 : g2); push!(δ2, isnan(d2) ? 0.0 : d2)
        push!(η, isnan(et) ? 0.0 : et); push!(νVC, isnan(nv) ? 0.0 : nv); push!(Y, y)
        push!(mm, molar_mass(m, i))
    end
    p = sortperm(ν0)
    return LineDatabase(; mol = mol[p], iso = iso[p], ν0 = FT.(ν0[p]), S = FT.(S[p]),
                        E_lower = FT.(E[p]), g_upper = FT.(gup[p]), γ_air = FT.(γair[p]),
                        n_self = FT.(nself[p]), δ_self = FT.(δself[p]),
                        γ_self = FT.(γself[p]), n_air = FT.(nair[p]), δ_air = FT.(δair[p]),
                        γ2_air = FT.(γ2[p]), δ2_air = FT.(δ2[p]), η = FT.(η[p]), νVC = FT.(νVC[p]),
                        Y_LM = FT.(Y[p]), molar_mass = FT.(mm[p]),
                        meta = SourceMetadata(source, T_REF, P_REF), partition = TIPS2021PF())
end
