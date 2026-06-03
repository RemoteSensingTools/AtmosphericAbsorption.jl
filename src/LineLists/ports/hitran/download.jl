#=
HITRAN direct download (HAPI's `fetch()` equivalent). Pulls standard 160-char `.par`
line data (Voigt parameters) from hitran.org's public line-by-line API into a Scratch
cache with `.meta.toml` provenance (edition, window, source URL, SHA-256, date). The
cache is reused only when it covers the requested wavenumber window.

Non-Voigt (Hartmann-Tran / speed-dependent / line-mixing) parameters require the
authenticated HITRANonline API (HAPI2) and a user API key — added separately.
=#

const _HITRAN_LBL_API = "https://hitran.org/lbl/api"

_hitran_dir() = @get_scratch!("hitran")

# Global isotopologue ids for a HITRAN molecule (from the iso_info table).
_hitran_iso_ids(mol::Integer) =
    [Int(_ISO_TABLE.global_id[k]) for k in eachindex(_ISO_TABLE.mol) if _ISO_TABLE.mol[k] == mol]

# Does the first non-empty line look like a 160-char HITRAN record (not an error page)?
function _looks_like_par(path)
    return open(path) do io
        for ln in eachline(io)
            isempty(strip(ln)) && continue
            return length(ln) ≥ 160 && tryparse(Int, strip(ln[1:2])) !== nothing
        end
        return false
    end
end

# Does the cached metadata cover [numin, numax]?
function _hitran_cache_covers(meta_path, numin, numax)
    isfile(meta_path) || return false
    meta = try
        TOML.parsefile(meta_path)
    catch
        return false
    end
    lo, hi = get(meta, "numin", nothing), get(meta, "numax", nothing)
    return lo isa Real && hi isa Real && lo ≤ numin && hi ≥ numax
end

"""
    fetch_hitran(molecule; numin=0, numax=150000, edition="HITRAN2020", force=false) -> path

Download HITRAN `.par` line data for `molecule` (all isotopologues) over `[numin, numax]`
cm⁻¹ from hitran.org, caching it in a Scratch dir with provenance. Returns the local
`.par` path. `edition` is the cache label; the public API serves the current HITRAN.
Wrap the result in a `HitranPort` to load it.
"""
function fetch_hitran(molecule::AbstractString; numin::Real = 0, numax::Real = 150000,
                      edition::AbstractString = "HITRAN2020", force::Bool = false)
    ids = _hitran_iso_ids(molecule_number(molecule))
    isempty(ids) && error("no HITRAN isotopologues found for \"$molecule\"")
    dir = joinpath(_hitran_dir(), edition)
    mkpath(dir)
    par, meta = joinpath(dir, "$molecule.par"), joinpath(dir, "$molecule.meta.toml")

    if !force && isfile(par) && _hitran_cache_covers(meta, numin, numax)
        return par
    end

    url = "$(_HITRAN_LBL_API)?iso_ids_list=$(join(ids, ','))&numin=$(numin)&numax=$(numax)"
    # hitran.org uses chunked transfer with no Content-Length, so Downloads reports a
    # spurious "bytes missing" error on success — request with throw=false, then validate
    # the HTTP status AND the content (so an error page is never cached as a .par).
    resp = open(par, "w") do out
        Downloads.request(url; output = out, throw = false)
    end
    status = resp isa Downloads.Response ? resp.status :
             (resp isa Downloads.RequestError && resp.response isa Downloads.Response) ?
             resp.response.status : 0
    if status == 403
        rm(par; force = true)
        error("HITRAN API rate limit exceeded (HTTP 403) for $molecule — try again later.")
    elseif status ≥ 400
        rm(par; force = true)
        error("HITRAN API request failed (HTTP $status) for $molecule.")
    end
    (isfile(par) && filesize(par) > 0 && _looks_like_par(par)) ||
        (rm(par; force = true);
         error("HITRAN download for $molecule [$numin, $numax] returned no valid line data."))

    open(meta, "w") do io
        println(io, "molecule = \"", molecule, "\"")
        println(io, "edition = \"", edition, "\"")
        println(io, "numin = ", numin)
        println(io, "numax = ", numax)
        println(io, "source = \"", url, "\"")
        println(io, "sha256 = \"", bytes2hex(SHA.sha256(read(par))), "\"")
        println(io, "fetched = \"", Dates.today(), "\"")
    end
    return par
end

"""
    HitranPort(; molecule, numin=0, numax=150000, edition="HITRAN2020", force=false)

Convenience constructor that downloads HITRAN data for `molecule` (via `fetch_hitran`)
and returns a `HitranPort` pointing at the cached `.par`.
"""
HitranPort(; molecule::AbstractString, numin::Real = 0, numax::Real = 150000,
           edition::AbstractString = "HITRAN2020", force::Bool = false) =
    HitranPort(fetch_hitran(molecule; numin, numax, edition, force); edition)
