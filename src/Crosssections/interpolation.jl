#=
Cross-section interpolation table: precompute a line-by-line model's σ(ν) on a (pressure,
temperature) grid once, then look it up by fast multilinear interpolation instead of re-summing
lines. The line-by-line cost is paid once at build time; queries are O(grid) and do no line work.
The table is architecture-aware: its σ cube (and ν nodes) live on `array_type(architecture)`, the
(p,T) bracket + corner blend is an arch-generic broadcast, and the ν-resample is a KernelAbstractions
kernel — so a query runs on CPU or GPU and returns an array on the model's architecture.
=#

"""
    InterpolationModel{FT}

A tabulated cross-section cube `σ[iν, ip, iT]` built from a `LineByLineModel`, with its node
vectors `ν` [cm⁻¹], `p` [hPa], `T` [K], the `vmr` it was built at, and the `architecture` it lives
on. `compute_cross_section(im, grid, p, T)` interpolates it (bilinear in (p, T), linear in ν, clamped
in p/T and zero outside the ν range) on that architecture — fast, no line-by-line work. The `σ`/`ν`
arrays are CPU `Array`s for `CPU()` and device arrays for `GPU()`/`MetalGPU()`. The table is
vmr-specific (`vmr` is fixed at build time). Build it with [`build_interpolation_model`](@ref);
persist with [`save_interpolation_model`](@ref).
"""
struct InterpolationModel{FT<:AbstractFloat,V<:AbstractVector{FT},
                          C<:AbstractArray{FT,3},A<:AbstractArchitecture} <: AbstractCrossSectionModel
    ν::V               # wavenumber nodes [cm⁻¹] (on the architecture)
    p::Vector{FT}      # pressure nodes [hPa], ascending (host — bracketed on host)
    T::Vector{FT}      # temperature nodes [K], ascending (host)
    σ::C               # σ[iν, ip, iT] [cm²/molecule] (on the architecture)
    vmr::FT            # broadener volume mixing ratio fixed at build time
    architecture::A
end

"""
    build_interpolation_model(model::LineByLineModel, ν, pressures, temperatures;
                              vmr=model.vmr, architecture=model.architecture) -> InterpolationModel

Evaluate the line-by-line `model` on wavenumber grid `ν` at every (pressure, temperature) node and
store the cube for fast interpolation later. `pressures` [hPa] and `temperatures` [K] are the ascending
interpolation nodes; the (expensive) line-by-line summation runs once, here. The resulting table lives
on `architecture` (defaults to the model's), so queries run there. Choose nodes fine enough that
interpolation between them is accurate for your use.
"""
function build_interpolation_model(model::LineByLineModel{FT}, ν::AbstractVector,
                                   pressures, temperatures; vmr::Real = model.vmr,
                                   architecture::AbstractArchitecture = model.architecture) where {FT}
    νv, ps, Ts = collect(FT, ν), collect(FT, pressures), collect(FT, temperatures)
    (issorted(ps) && issorted(Ts)) ||
        throw(ArgumentError("pressures and temperatures must be ascending"))
    σh = Array{FT,3}(undef, length(νv), length(ps), length(Ts))   # assemble on host…
    for (k, T) in enumerate(Ts), (j, p) in enumerate(ps)
        σh[:, j, k] .= Array(compute_cross_section(model, νv, p, T; vmr))
    end
    AT = array_type(architecture)                                  # …then place on the architecture
    return InterpolationModel(AT(νv), ps, Ts, AT(σh), FT(vmr), architecture)
end

# Lower bracket index + interpolation fraction of x in ascending `nodes`, clamped so that both `i`
# and `i+1` are valid (single-node grids return fraction 0). Host-side scalar helper (p/T/vmr axes).
@inline function _bracket(nodes, x)
    n = length(nodes)
    n == 1 && return 1, 1, zero(eltype(nodes))
    i = clamp(searchsortedlast(nodes, x), 1, n - 1)
    return i, i + 1, (x - nodes[i]) / (nodes[i+1] - nodes[i])
end

# Zero outside the tabulated ν range (matching TabulatedCrossSection), linear within — host helper.
@inline function _interp_linear(xs, ys, x)
    (x < first(xs) || x > last(xs)) && return zero(eltype(ys))
    i, i1, f = _bracket(xs, x)
    return (1 - f) * ys[i] + f * ys[i1]
end

# Resample σ(ν) given on the table's `ν` nodes onto `grid`: linear between nodes, zero outside the
# range. One KernelAbstractions kernel serves CPU and GPU; each work-item does one output point
# (binary search in `ν` + linear blend). Shared by InterpolationModel and AbscoLUT.
@kernel function _lut_resample_kernel!(out, @Const(grid), @Const(ν), @Const(σν))
    i = @index(Global)
    @inbounds begin
        FT = eltype(out)
        x = grid[i]
        n = length(ν)
        if x < ν[1] || x > ν[n]
            out[i] = zero(FT)
        elseif n == 1
            out[i] = σν[1]
        else
            lo = 1
            hi = n
            while hi - lo > 1                       # last index with ν[lo] ≤ x
                mid = (lo + hi) >>> 1
                ν[mid] <= x ? (lo = mid) : (hi = mid)
            end
            f = (x - ν[lo]) / (ν[lo+1] - ν[lo])
            out[i] = (1 - f) * σν[lo] + f * σν[lo+1]
        end
    end
end

# σν is already on `arch`; resample onto `grid` (also moved to `arch`) and return an `arch` array.
function _lut_resample(arch::AbstractArchitecture, ν, σν, grid)
    FT = eltype(σν)
    Ng = length(grid)
    out = array_type(arch)(Vector{FT}(undef, Ng))
    Ng == 0 && return out
    gridd = array_type(arch)(collect(FT, grid))
    _lut_resample_kernel!(devi(arch))(out, gridd, ν, σν; ndrange = Ng)
    synchronize_if_gpu(arch)
    return out
end

"""
    compute_cross_section(im::InterpolationModel, grid, pressure, temperature) -> array

Cross-section [cm²/molecule] on `grid` [cm⁻¹] at `pressure` [hPa], `temperature` [K], bilinearly
interpolated in (p, T) — clamped to the tabulated range — then linearly resampled in ν (zero outside
the table's ν range). `vmr` is fixed at build time. Runs on the table's architecture and returns an
array there (host `Array` for `CPU()`, device array for a GPU). Querying on the table's own `ν` skips
the ν resample.
"""
function compute_cross_section(im::InterpolationModel{FT}, grid::AbstractVector,
                               pressure::Real, temperature::Real) where {FT}
    jp, jp1, fp = _bracket(im.p, clamp(FT(pressure), first(im.p), last(im.p)))
    kt, kt1, ft = _bracket(im.T, clamp(FT(temperature), first(im.T), last(im.T)))
    # σ(ν) on the table grid by an arch-generic bilinear blend over the four (p, T) corner panels.
    σν = @views (1 - fp) * (1 - ft) .* im.σ[:, jp, kt]  .+ fp * (1 - ft) .* im.σ[:, jp1, kt] .+
                (1 - fp) * ft       .* im.σ[:, jp, kt1] .+ fp * ft       .* im.σ[:, jp1, kt1]
    grid === im.ν && return σν                            # fast path: query is the table grid
    return _lut_resample(im.architecture, im.ν, σν, grid)
end

"""
    save_interpolation_model(path, im)
    load_interpolation_model(path) -> InterpolationModel

Persist / restore an interpolation table to/from `path` (via Julia `Serialization`). Saving copies
the table to the host (`CPU()`), so the file is portable; move it back with `architecture=` on a
rebuild if you need it on a GPU. The format is tied to the Julia version that wrote it — rebuild the
table after a Julia upgrade.
"""
function save_interpolation_model(path::AbstractString, im::InterpolationModel)
    cpu = InterpolationModel(Array(im.ν), im.p, im.T, Array(im.σ), im.vmr, Architectures.CPU())
    open(io -> Serialization.serialize(io, cpu), path, "w")
end
load_interpolation_model(path::AbstractString) = open(Serialization.deserialize, path)
