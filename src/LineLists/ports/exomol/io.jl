#=
ExoMol file fetching and parsing. Files live at
www.exomol.com/db/<molecule>/<iso_slug>/<linelist>/<iso_slug>__<linelist>.<ext>
with the shared broadening one level up (<iso_slug>__air.broad). `.states`/`.trans`
are bzip2-compressed and streamed; `.def`/`.pf`/`.broad` are plain text. Large line
lists split `.trans` into wavenumber chunks (`__NNNNN-NNNNN.trans.bz2`, bin =
max_wavenumber / n_files) so a windowed query downloads only the chunks it touches.

Downloads land in a Scratch cache with a `.meta.toml` recording the dataset version
(from the `.def`), source URL, `.def` SHA-256, and fetch date, so runs are traceable.
All formats are whitespace-delimited (field widths vary; fixed slicing is unsafe).
Verified against CO/12C-16O/Li2015.
=#

const _EXOMOL_BASE = "https://www.exomol.com/db"

_exomol_dir() = @get_scratch!("exomol")
_base(port)   = "$(_EXOMOL_BASE)/$(port.molecule)/$(port.iso_slug)/$(port.linelist)"
_isobase(port)= "$(_EXOMOL_BASE)/$(port.molecule)/$(port.iso_slug)"
_stem(port)   = "$(port.iso_slug)__$(port.linelist)"

# Stream `path`, transparently bzip2-decompressing `.bz2`, and call `f(io)`.
function _open_text(f, path)
    if endswith(path, ".bz2")
        io = Bzip2DecompressorStream(open(path))
        try
            return f(io)
        finally
            close(io)
        end
    else
        return open(f, path)
    end
end

function _cached(dir, name, url; optional = false)
    dest = joinpath(dir, name)
    isfile(dest) && return dest
    try
        Downloads.download(url, dest)
    catch e
        optional && return dest          # e.g. a linelist with no .broad
        rethrow(e)
    end
    return dest
end

# Record provenance (version, source, .def hash, date) next to the cached files.
function _write_provenance(dir, port, def_path, version)
    meta = joinpath(dir, "$(_stem(port)).meta.toml")
    isfile(meta) && return meta
    sha = bytes2hex(SHA.sha256(read(def_path)))
    open(meta, "w") do io
        println(io, "version = \"", version, "\"")
        println(io, "source = \"", _base(port), "\"")
        println(io, "def_sha256 = \"", sha, "\"")
        println(io, "fetched = \"", Dates.today(), "\"")
    end
    return meta
end

"""
    fetch_exomol_meta(port) -> NamedTuple

Download (and cache) the small ExoMol files — `.def`, `.states.bz2`, `.pf`, and the
isotopologue `.broad` — and write provenance. The (possibly chunked) `.trans` files
are fetched separately by `fetch_trans` once the wavenumber window is known.
"""
function fetch_exomol_meta(port)
    dir  = _exomol_dir()
    base, isob, stem = _base(port), _isobase(port), _stem(port)
    def    = _cached(dir, "$stem.def",        "$base/$stem.def")
    states = _cached(dir, "$stem.states.bz2", "$base/$stem.states.bz2")
    pf     = _cached(dir, "$stem.pf",         "$base/$stem.pf")
    broad  = _cached(dir, "$(port.iso_slug)__air.broad",
                     "$isob/$(port.iso_slug)__air.broad"; optional = true)
    return (; dir, def, states, pf, broad)
end

"""
    fetch_pf(port) -> path

Download only the `.pf` partition-function file (avoids pulling multi-GB transitions
when a caller just needs Q(T)).
"""
fetch_pf(port) = _cached(_exomol_dir(), "$(_stem(port)).pf", "$(_base(port))/$(_stem(port)).pf")

"""
    fetch_trans(port, def, ν_min, ν_max) -> Vector{String}

Download and return the local paths of the `.trans` chunk(s) covering `[ν_min, ν_max]`.
Single-file line lists return one path; chunked lists return only the overlapping bins.
"""
function fetch_trans(port, def, ν_min, ν_max)
    dir, base, stem = _exomol_dir(), _base(port), _stem(port)
    if def.ntransfiles ≤ 1
        return [_cached(dir, "$stem.trans.bz2", "$base/$stem.trans.bz2")]
    end
    bin  = round(Int, def.maxwav / def.ntransfiles)
    lo   = max(0, floor(Int, ν_min / bin) * bin)
    hi   = min(def.ntransfiles * bin, isfinite(ν_max) ? ceil(Int, ν_max / bin) * bin : def.ntransfiles * bin)
    paths = String[]
    for s in lo:bin:(hi - 1)
        name = @sprintf("%s__%05d-%05d.trans.bz2", stem, s, s + bin)
        push!(paths, _cached(dir, name, "$base/$name"))
    end
    return paths
end

"""
    parse_def(path) -> NamedTuple

Read the fields the cross-section needs from an ExoMol `.def` (matched by comment
text, robust to the conditional quantum-label block): dataset `version`, isotopologue
molar `mass` [g/mol], state count `nstates`, transition-file count `ntransfiles` and
maximum wavenumber `maxwav` (for chunk windowing), and the default broadening.
"""
function parse_def(path)
    version = "unknown"
    mass = 0.0; nstates = 0; ntransfiles = 1; maxwav = 0.0
    default_γ = 0.07; default_n = 0.5
    for ln in eachline(path)
        parts = split(ln, '#'; limit = 2)
        length(parts) < 2 && continue
        val, desc = strip(parts[1]), lowercase(parts[2])
        isempty(val) && continue
        if occursin("version number", desc)
            version = String(first(split(val)))
        elseif occursin("isotopologue mass", desc)
            mass = parse(Float64, first(split(val)))          # "Da kg" → take Da (g/mol)
        elseif occursin("no. of states", desc)
            nstates = parse(Int, first(split(val)))
        elseif occursin("no. of transition files", desc)
            ntransfiles = parse(Int, first(split(val)))
        elseif occursin("maximum wavenumber", desc)
            maxwav = parse(Float64, first(split(val)))
        elseif occursin("default value of lorentzian half-width", desc)
            default_γ = parse(Float64, first(split(val)))
        elseif occursin("default value of temperature exponent", desc)
            default_n = parse(Float64, first(split(val)))
        end
    end
    return (; version, mass, nstates, ntransfiles, maxwav, default_γ, default_n)
end

"""
    read_states(path, nstates) -> NamedTuple

Load an ExoMol `.states` file into dense, 1-indexed arrays `(E[cm⁻¹], g, J)` keyed by
the contiguous state id (so `E[id]` is O(1)). `J` is `Float64` to accept half-integer
rotational quanta (open-shell molecules).
"""
function read_states(path, nstates)
    E = Vector{Float64}(undef, nstates)
    g = Vector{Float64}(undef, nstates)
    J = Vector{Float64}(undef, nstates)
    _open_text(path) do io
        for ln in eachline(io)
            t = split(ln)
            isempty(t) && continue
            i = parse(Int, t[1])
            1 ≤ i ≤ nstates || error("ExoMol .states id $i outside 1:$nstates ($path)")
            E[i] = parse(Float64, t[2]); g[i] = parse(Float64, t[3]); J[i] = parse(Float64, t[4])
        end
    end
    return (; E, g, J)
end

"""
    read_pf(path) -> (T, Q)

Read an ExoMol `.pf` partition-function table into temperature/partition vectors.
"""
function read_pf(path)
    T = Float64[]; Q = Float64[]
    for ln in eachline(path)
        t = split(ln)
        isempty(t) && continue
        push!(T, parse(Float64, t[1])); push!(Q, parse(Float64, t[2]))
    end
    return (; T, Q)
end

"""
    read_broad(path) -> (diet, table)

Parse an ExoMol `.broad` file into a `key → (γ_ref[cm⁻¹/bar], n)` lookup. `diet` is the
first-column code: `"a0"` keys on lower-state J, `"m0"` on the signed branch index m.
Returns `("", empty)` if the file is absent (caller falls back to the `.def` defaults).
"""
function read_broad(path)
    table = Dict{Int,Tuple{Float64,Float64}}()
    diet  = ""
    isfile(path) || return (diet, table)
    for ln in eachline(path)
        t = split(ln)
        length(t) < 4 && continue
        diet = t[1]
        table[parse(Int, t[4])] = (parse(Float64, t[2]), parse(Float64, t[3]))
    end
    return (diet, table)
end
