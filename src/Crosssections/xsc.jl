#=
HITRAN tabulated absorption cross-sections (.xsc) — measured σ(ν) panels at fixed
(T, p) for heavy / IR-active molecules (CFCs, SF₆, VOCs, …) that have no usable line
list. Parsed from HITRAN's semi-fixed-width format (a 100-char header followed by
`npnts` σ values, ten per line, each ten characters wide), one or more (T, p) panels
per file, into `XscBand`s. `TabulatedCrossSection` interpolates them onto a requested
grid: linear in wavenumber, linear in temperature at the nearest tabulated pressure —
the same read→interpolate→project shape as the continuum readers.
=#

const _TORR_TO_HPA = 101325 / 760 / 100   # 1 Torr in hPa

"""
    XscBand{FT}

One tabulated cross-section panel: σ(ν) sampled at `length(σ)` points evenly spanning
[`νmin`, `νmax`] cm⁻¹, measured at temperature `T` [K] and pressure `p` [hPa].
"""
struct XscBand{FT<:AbstractFloat}
    molecule::String
    T::FT          # K
    p::FT          # hPa (converted from the file's Torr)
    νmin::FT
    νmax::FT
    σ::Vector{FT}  # cm²/molecule
end

# Linear interpolation of a panel's σ at wavenumber ν; 0 outside [νmin, νmax].
@inline function _band_sigma(b::XscBand{FT}, ν::FT) where {FT}
    (ν < b.νmin || ν > b.νmax) && return zero(FT)
    last = length(b.σ) - 1
    x = clamp((ν - b.νmin) / (b.νmax - b.νmin) * last, zero(FT), FT(last))  # fractional 0-based index
    i = unsafe_trunc(Int, x)
    i ≥ last && return b.σ[end]
    f = x - i
    return (1 - f) * b.σ[i + 1] + f * b.σ[i + 2]              # σ is 1-based
end

"""
    read_xsc(path; FT=Float64) -> Vector{XscBand{FT}}

Parse a HITRAN cross-section (`.xsc`) file — one or more (T, p) panels — into `XscBand`s.
"""
function read_xsc(path::AbstractString; FT::Type{<:AbstractFloat} = Float64)
    bands = XscBand{FT}[]
    open(path) do io
        while !eof(io)
            header = readline(io)
            length(header) < 60 && continue                   # blank/short line between panels
            mol   = String(strip(header[1:20]))
            νmin  = parse(Float64, strip(header[21:30]))
            νmax  = parse(Float64, strip(header[31:40]))
            npnts = parse(Int, strip(header[41:47]))
            T     = parse(Float64, strip(header[48:54]))
            p     = parse(Float64, strip(header[55:60])) * _TORR_TO_HPA
            σ = Vector{FT}(undef, npnts)
            k = 0
            while k < npnts                                   # ten fixed 10-char fields per line
                line = readline(io)
                c = 1
                while c ≤ length(line) && k < npnts
                    tok = strip(SubString(line, c, min(c + 9, length(line))))
                    isempty(tok) || (σ[k += 1] = parse(FT, tok))
                    c += 10
                end
            end
            push!(bands, XscBand{FT}(mol, FT(T), FT(p), FT(νmin), FT(νmax), σ))
        end
    end
    return bands
end

"""
    TabulatedCrossSection{FT}

A pre-computed cross-section source for one molecule, holding the `.xsc` panels over
(wavenumber, temperature, pressure). `compute_cross_section(model, grid, p, T)`
interpolates σ onto `grid` — linear in ν, linear in T at the nearest tabulated p (full
(T, p) bilinear interpolation is not done). It always returns a CPU `Vector{FT}`; wrap
it with `array_type(arch)(...)` to move it onto a device.
"""
struct TabulatedCrossSection{FT<:AbstractFloat} <: AbstractCrossSectionModel
    molecule::String
    bands::Vector{XscBand{FT}}
end

"""
    load_xsc(path; FT=Float64) -> TabulatedCrossSection

Read a HITRAN `.xsc` file into a tabulated cross-section model.
"""
function load_xsc(path::AbstractString; FT::Type{<:AbstractFloat} = Float64)
    bands = read_xsc(path; FT)
    isempty(bands) && error("no cross-section panels parsed from $path")
    return TabulatedCrossSection{FT}(bands[1].molecule, bands)
end

# σ(ν; p, T): among the panels covering ν, pick the nearest tabulated pressure and
# interpolate linearly in temperature (clamped to the tabulated range). Returns 0 where
# no panel covers ν.
function _interp_xsc(bands::Vector{XscBand{FT}}, ν::FT, p::FT, T::FT) where {FT}
    pbest = zero(FT); dpmin = FT(Inf); covered = false
    for b in bands
        b.νmin ≤ ν ≤ b.νmax || continue
        covered = true
        if abs(b.p - p) < dpmin
            dpmin = abs(b.p - p); pbest = b.p
        end
    end
    covered || return zero(FT)
    # Bracket T among the covering panels at the chosen pressure. At least one such panel
    # exists (covered), so haslo || hashi always holds; a single bracket is T-extrapolation.
    Tlo = Thi = T; σlo = σhi = zero(FT); haslo = hashi = false
    for b in bands
        (b.p == pbest && b.νmin ≤ ν ≤ b.νmax) || continue
        s = _band_sigma(b, ν)
        if b.T ≤ T && (!haslo || b.T > Tlo)
            Tlo, σlo, haslo = b.T, s, true
        end
        if b.T ≥ T && (!hashi || b.T < Thi)
            Thi, σhi, hashi = b.T, s, true
        end
    end
    haslo && hashi || return haslo ? σlo : σhi
    Thi == Tlo && return σlo
    return σlo + (σhi - σlo) * (T - Tlo) / (Thi - Tlo)
end

"""
    compute_cross_section(model::TabulatedCrossSection, grid, p, T) -> Vector

Cross-section [cm²/molecule] on `grid` [cm⁻¹] at pressure `p` [hPa], temperature `T` [K],
interpolated from the tabulated `.xsc` panels.
"""
function compute_cross_section(model::TabulatedCrossSection{FT}, grid::AbstractVector,
                               pressure::Real, temperature::Real) where {FT}
    p, T = FT(pressure), FT(temperature)
    return FT[_interp_xsc(model.bands, FT(ν), p, T) for ν in grid]
end

"""
    fetch_hitran_xsc(filename; force=false) -> path

Download a HITRAN cross-section file by its `.xsc` filename (as listed on hitran.org)
into the scratch cache, returning the local path. Reused on subsequent calls.
"""
function fetch_hitran_xsc(filename::AbstractString; force::Bool = false)
    dir = @get_scratch!("hitran_xsec"); path = joinpath(dir, filename)
    (force || !isfile(path) || filesize(path) == 0) || return path
    Downloads.download("https://hitran.org/data/xsec/$filename", path)
    (isfile(path) && filesize(path) > 0) ||
        (rm(path; force = true); error("HITRAN xsc download failed for $filename"))
    return path
end
