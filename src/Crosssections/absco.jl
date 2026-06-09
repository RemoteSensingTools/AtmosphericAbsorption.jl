#=
ABSCO lookup table. AER's ABSCO files tabulate σ(ν, p, T, H₂O-broadener) on a grid whose temperature
nodes *slide with pressure* (T is stored per-pressure) and whose H₂O broadener VMR is a third axis.
Rather than resample onto a regular (p, T) cube — which wastes storage on the unphysical
(low-p/high-T) corners ABSCO never covers and loses exact recoverability — we keep the native grid
and interpolate at query time: bracket p, interpolate in each bracketing pressure's own T-axis, blend
in p and in the broadener VMR, then linear in ν. Querying at an original node returns the stored value
exactly. Read real ABSCO files with `read_absco` (needs NCDatasets). Like the other cross-section
models it is architecture-aware: the σ cube / ν nodes live on `array_type(architecture)`, the blend is
an arch-generic broadcast and the ν-resample is the shared KernelAbstractions kernel — so a query runs
on CPU or GPU and returns an array there.
=#

"""
    AbscoLUT{FT}

AER ABSCO tabulated cross-sections on their native grid: `σ[iν, ivmr, iT, ip]` [cm²/molecule] with
wavenumber nodes `ν` [cm⁻¹], ascending pressure nodes `p` [hPa], **per-pressure** temperature nodes
`T[iT, ip]` [K] (the ABSCO temperature grid slides with pressure), and H₂O broadener VMR nodes `vmr`.
Query with `compute_cross_section(lut, grid, p, T; vmr, interp)`, interpolated in (p, T, vmr) and
linear in ν — querying at an original node is exact. The σ/ν arrays live on `architecture`. Build from
a file with [`read_absco`](@ref).
"""
struct AbscoLUT{FT<:AbstractFloat,V<:AbstractVector{FT},
                C<:AbstractArray{FT,4},A<:AbstractArchitecture} <: AbstractCrossSectionModel
    mol::Int
    iso::Int
    ν::V               # wavenumber nodes [cm⁻¹] (on the architecture)
    p::Vector{FT}      # pressure nodes [hPa], ascending (host — bracketed on host)
    T::Matrix{FT}      # T[iT, ip] temperature nodes [K] per pressure, ascending in iT (host)
    vmr::Vector{FT}    # H₂O broadener VMR nodes, ascending (host)
    σ::C               # σ[iν, ivmr, iT, ip] [cm²/molecule] (on the architecture)
    architecture::A
end

Base.eltype(::AbscoLUT{FT}) where {FT} = FT

"""
    AbscoLUT(mol, iso, ν, p, T, vmr, σ; architecture=default_architecture())

Assemble an [`AbscoLUT`](@ref), placing the wavenumber nodes `ν` and the cube `σ[iν, ivmr, iT, ip]`
on `architecture` (the per-pressure `T[iT, ip]`, `p` and broadener `vmr` axes stay on the host, where
they are bracketed). `read_absco` calls this after reading the file.
"""
function AbscoLUT(mol::Integer, iso::Integer, ν::AbstractVector{FT}, p::AbstractVector{FT},
                  T::AbstractMatrix{FT}, vmr::AbstractVector{FT}, σ::AbstractArray{FT,4};
                  architecture::AbstractArchitecture = default_architecture()) where {FT}
    AT = array_type(architecture)
    return AbscoLUT(Int(mol), Int(iso), AT(Vector{FT}(ν)), Vector{FT}(p), Matrix{FT}(T),
                    Vector{FT}(vmr), AT(σ), architecture)
end

# σ(ν) at (pressure node ip, broadener node kvi), interpolated in temperature between index it and
# it1 (fraction ft). `:cubic` uses a uniform 4-point Catmull-Rom in the T-index (ABSCO's per-pressure
# T nodes are equally spaced) — node-exact like linear but smooth; `:linear` is the two-node blend.
# Pure broadcast over (device) views, so it runs on CPU or GPU.
@inline function _absco_interp_T(σ::AbstractArray{FT,4}, ip::Int, kvi::Int, it::Int, it1::Int,
                                 ft::FT, interp::Symbol) where {FT}
    if interp === :cubic
        nT = size(σ, 3)
        y0 = @view σ[:, kvi, max(it - 1, 1), ip]
        y1 = @view σ[:, kvi, it, ip]
        y2 = @view σ[:, kvi, it1, ip]
        y3 = @view σ[:, kvi, min(it1 + 1, nT), ip]
        return @. y1 + ft * (FT(0.5) * (y2 - y0) +
                  ft * ((y0 - FT(2.5) * y1 + FT(2) * y2 - FT(0.5) * y3) +
                  ft * (FT(0.5) * (y3 - y0) + FT(1.5) * (y1 - y2))))
    else
        y1 = @view σ[:, kvi, it, ip]
        y2 = @view σ[:, kvi, it1, ip]
        return @. (1 - ft) * y1 + ft * y2
    end
end

# σ(ν) at pressure node ip and temperature T, blended over the broadener bracket (kv, kv1, fv).
function _absco_at_p(lut::AbscoLUT{FT}, ip::Int, T::FT, kv::Int, kv1::Int, fv::FT,
                     interp::Symbol) where {FT}
    Taxis = @view lut.T[:, ip]
    it, it1, ft = _bracket(Taxis, clamp(T, first(Taxis), last(Taxis)))
    σ0 = _absco_interp_T(lut.σ, ip, kv, it, it1, ft, interp)
    kv == kv1 && return σ0
    σ1 = _absco_interp_T(lut.σ, ip, kv1, it, it1, ft, interp)
    return @. (1 - fv) * σ0 + fv * σ1
end

"""
    compute_cross_section(lut::AbscoLUT, grid, pressure, temperature; vmr=0, interp=:linear) -> array

Cross-section [cm²/molecule] on `grid` [cm⁻¹] at `pressure` [hPa], `temperature` [K] and H₂O broadener
`vmr`, looked up from the ABSCO table: interpolated in pressure, in each bracketing pressure's own
temperature axis, and in the broadener VMR (all clamped to the tabulated range), then linear in ν (zero
outside the table's ν range). `interp` selects the **temperature** interpolation — `:linear` (default)
or `:cubic` (smooth Catmull-Rom); both are exact at nodes. Runs on the table's architecture and returns
an array there; querying on the table's own `ν` skips the ν resample.
"""
function compute_cross_section(lut::AbscoLUT{FT}, grid::AbstractVector, pressure::Real,
                               temperature::Real; vmr::Real = 0, interp::Symbol = :linear) where {FT}
    jp, jp1, fp = _bracket(lut.p, clamp(FT(pressure), first(lut.p), last(lut.p)))
    kv, kv1, fv = _bracket(lut.vmr, clamp(FT(vmr), first(lut.vmr), last(lut.vmr)))
    T = FT(temperature)
    σlo = _absco_at_p(lut, jp, T, kv, kv1, fv, interp)
    if jp == jp1
        σν = σlo
    else
        σhi = _absco_at_p(lut, jp1, T, kv, kv1, fv, interp)
        σν = @. (1 - fp) * σlo + fp * σhi
    end
    grid === lut.ν && return σν
    return _lut_resample(lut.architecture, lut.ν, σν, grid)
end

"""
    read_absco(path; scale=1.0, FT=Float32, architecture=default_architecture()) -> AbscoLUT

Read an AER ABSCO `.hdf` file into an [`AbscoLUT`](@ref) on `architecture`. **Requires NCDatasets**
(`using NCDatasets` to enable it — the `.hdf` files are HDF5, which NCDatasets reads natively); the
method lives in a package extension so AtmosphericAbsorption's core stays NetCDF-free.
"""
function read_absco end

"""
    save_absco_lut(path, lut)
    load_absco_lut(path) -> AbscoLUT

Persist / restore an [`AbscoLUT`](@ref) via Julia `Serialization` (copied to the host on save, so the
file is portable; format tied to the writing Julia version). Use this to ship the converted table
without the raw ABSCO file.
"""
function save_absco_lut(path::AbstractString, lut::AbscoLUT)
    cpu = AbscoLUT(lut.mol, lut.iso, Array(lut.ν), lut.p, lut.T, lut.vmr, Array(lut.σ),
                   Architectures.CPU())
    open(io -> Serialization.serialize(io, cpu), path, "w")
end
load_absco_lut(path::AbstractString) = open(Serialization.deserialize, path)
