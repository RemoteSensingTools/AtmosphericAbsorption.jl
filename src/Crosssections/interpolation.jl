#=
Cross-section interpolation table: precompute a line-by-line model's σ(ν) on a (pressure,
temperature) grid once, then look it up by fast multilinear interpolation instead of
re-summing lines. The line-by-line cost is paid once at build time; queries are O(grid) and
do no line work — the same precomputed-LUT idea the legacy Absorption module used to feed RT.
=#

"""
    InterpolationModel{FT}

A tabulated cross-section cube `σ[iν, ip, iT]` built from a `LineByLineModel`, with its node
vectors `ν` [cm⁻¹], `p` [hPa], `T` [K] and the `vmr` it was built at. `compute_cross_section(im,
grid, p, T)` interpolates it (bilinear in (p, T), linear in ν, clamped in p/T and zero outside
the ν range) — fast, no line-by-line work. The table is vmr-specific (`vmr` is fixed at build
time, not a per-call argument). Build it with [`build_interpolation_model`](@ref); persist with
[`save_interpolation_model`](@ref).
"""
struct InterpolationModel{FT<:AbstractFloat} <: AbstractCrossSectionModel
    ν::Vector{FT}      # wavenumber nodes [cm⁻¹]
    p::Vector{FT}      # pressure nodes [hPa], ascending
    T::Vector{FT}      # temperature nodes [K], ascending
    σ::Array{FT,3}     # σ[iν, ip, iT] [cm²/molecule]
    vmr::FT            # broadener volume mixing ratio fixed at build time
end

"""
    build_interpolation_model(model::LineByLineModel, ν, pressures, temperatures;
                              vmr=model.vmr) -> InterpolationModel

Evaluate the line-by-line `model` on wavenumber grid `ν` at every (pressure, temperature)
node and store the result for fast interpolation later. `pressures` [hPa] and `temperatures`
[K] are the ascending interpolation nodes; the (expensive) line-by-line summation runs once,
here. Choose the nodes fine enough that interpolation between them is accurate for your use.
"""
function build_interpolation_model(model::LineByLineModel{FT}, ν::AbstractVector,
                                   pressures, temperatures; vmr::Real = model.vmr) where {FT}
    νv, ps, Ts = collect(FT, ν), collect(FT, pressures), collect(FT, temperatures)
    (issorted(ps) && issorted(Ts)) ||
        throw(ArgumentError("pressures and temperatures must be ascending"))
    σ = Array{FT,3}(undef, length(νv), length(ps), length(Ts))
    for (k, T) in enumerate(Ts), (j, p) in enumerate(ps)
        σ[:, j, k] .= Array(compute_cross_section(model, νv, p, T; vmr))
    end
    return InterpolationModel(νv, ps, Ts, σ, FT(vmr))
end

# Lower bracket index + interpolation fraction of x in ascending `nodes`, clamped so that
# both `i` and `i+1` are valid (single-node grids return fraction 0).
@inline function _bracket(nodes, x)
    n = length(nodes)
    n == 1 && return 1, 1, zero(eltype(nodes))
    i = clamp(searchsortedlast(nodes, x), 1, n - 1)
    return i, i + 1, (x - nodes[i]) / (nodes[i+1] - nodes[i])
end

# Zero outside the tabulated ν range (no absorption recorded there), matching
# TabulatedCrossSection; inside, linear between the bracketing nodes.
@inline function _interp_linear(xs, ys, x)
    (x < first(xs) || x > last(xs)) && return zero(eltype(ys))
    i, i1, f = _bracket(xs, x)
    return (1 - f) * ys[i] + f * ys[i1]
end

"""
    compute_cross_section(im::InterpolationModel, grid, pressure, temperature) -> Vector

Cross-section [cm²/molecule] on `grid` [cm⁻¹] at `pressure` [hPa], `temperature` [K],
bilinearly interpolated in (p, T) — clamped to the tabulated p/T range — then linearly resampled
in ν (zero outside the table's ν range). `vmr` is fixed at build time, so it is not a per-call
argument here. Returns a CPU `Vector{FT}`; querying on the table's own `ν` is exact and allocation-light.
"""
function compute_cross_section(im::InterpolationModel{FT}, grid::AbstractVector,
                               pressure::Real, temperature::Real) where {FT}
    jp, jp1, fp = _bracket(im.p, clamp(FT(pressure), first(im.p), last(im.p)))
    kt, kt1, ft = _bracket(im.T, clamp(FT(temperature), first(im.T), last(im.T)))
    # σ(ν) on the table grid by bilinear blend over the four (p, T) corner panels.
    σν = @views (1 - fp) * (1 - ft) .* im.σ[:, jp, kt]  .+ fp * (1 - ft) .* im.σ[:, jp1, kt] .+
                (1 - fp) * ft       .* im.σ[:, jp, kt1] .+ fp * ft       .* im.σ[:, jp1, kt1]
    grid === im.ν && return collect(FT, σν)               # fast path: no ν resampling needed
    return FT[_interp_linear(im.ν, σν, FT(x)) for x in grid]
end

"""
    save_interpolation_model(path, im)
    load_interpolation_model(path) -> InterpolationModel

Persist / restore an interpolation table to/from `path` (via Julia `Serialization`). The
format is tied to the Julia version that wrote it — rebuild the table after a Julia upgrade.
"""
save_interpolation_model(path::AbstractString, im::InterpolationModel) =
    open(io -> Serialization.serialize(io, im), path, "w")
load_interpolation_model(path::AbstractString) = open(Serialization.deserialize, path)
