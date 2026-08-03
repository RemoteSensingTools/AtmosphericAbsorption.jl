"""Supertype for strategies that map a resolved spectrum to a requested grid."""
abstract type AbstractSpectralSampling end

"""
    PointSampling()

Evaluate a cross-section at the requested grid points. This is the default and
preserves the behavior of the original `compute_cross_section` API.
"""
struct PointSampling <: AbstractSpectralSampling end

abstract type AbstractConservativeSampling <: AbstractSpectralSampling end

struct ConservativeGridConfiguration{S,E}
    fine_step::S
    refinement::Int
    cell_edges::E
end

function ConservativeGridConfiguration(
    ;
    fine_step=nothing,
    refinement::Integer=32,
    cell_edges=nothing,
)
    refinement > 0 || throw(ArgumentError("refinement must be positive"))
    if fine_step !== nothing
        isfinite(fine_step) && fine_step > 0 || throw(ArgumentError(
            "fine_step must be finite and positive",
        ))
    end
    if cell_edges !== nothing
        all(isfinite, cell_edges) || throw(ArgumentError(
            "cell_edges must be finite",
        ))
    end
    return ConservativeGridConfiguration(
        fine_step,
        Int(refinement),
        cell_edges,
    )
end

"""
    ConservativeCrossSectionSampling(; fine_step=nothing, refinement=32,
                                       cell_edges=nothing)

Resolve a cross-section on a fine grid and return target-cell averages of
`σ`. The integral `sum(σ̄ .* diff(cell_edges))` is conserved. If `fine_step`
is omitted, the smallest target-cell width divided by `refinement` determines
the fine-grid spacing. Target-cell edges are inferred from grid centers unless
provided explicitly.
"""
struct ConservativeCrossSectionSampling{C<:ConservativeGridConfiguration} <:
       AbstractConservativeSampling
    configuration::C
end

function ConservativeCrossSectionSampling(; kwargs...)
    return ConservativeCrossSectionSampling(
        ConservativeGridConfiguration(; kwargs...),
    )
end

"""
    ConservativeTransmissionSampling(column_density_molecules_cm2;
                                      fine_step=nothing, refinement=32,
                                      cell_edges=nothing)

Resolve `σ` on a fine grid, conservatively average
`exp(-column_density_molecules_cm2 * σ)` in every target cell, and return the
effective cross-section `-log(T̄) / column_density_molecules_cm2`. Applying
Beer--Lambert at the same reference column therefore exactly reproduces the
cell-mean fine-grid transmission.

The result is column-dependent. Use [`ConservativeCrossSectionSampling`](@ref)
when a single cross-section must serve a wide range of columns, or conservatively
bin total transmission after combining all atmospheric layers.
"""
struct ConservativeTransmissionSampling{
    FT<:AbstractFloat,
    C<:ConservativeGridConfiguration,
} <: AbstractConservativeSampling
    column_density_molecules_cm2::FT
    configuration::C
end

function ConservativeTransmissionSampling(
    column_density_molecules_cm2::Real;
    kwargs...,
)
    column = float(column_density_molecules_cm2)
    isfinite(column) && column > 0 || throw(ArgumentError(
        "column_density_molecules_cm2 must be finite and positive",
    ))
    return ConservativeTransmissionSampling(
        column,
        ConservativeGridConfiguration(; kwargs...),
    )
end

sampling_configuration(sampling::AbstractConservativeSampling) =
    sampling.configuration

function _monotonic_direction(grid, label)
    length(grid) >= 2 || throw(ArgumentError(
        "$label requires at least two points",
    ))
    differences = diff(grid)
    all(differences .> zero(eltype(differences))) && return 1
    all(differences .< zero(eltype(differences))) && return -1
    throw(ArgumentError("$label must be strictly monotonic"))
end

"""
    spectral_cell_edges(grid)

Infer cell edges from strictly monotonic cell centers using midpoint interior
edges and half-spacing extrapolation at both endpoints. The returned edges have
the same ordering as `grid`.
"""
function spectral_cell_edges(grid::AbstractVector)
    _monotonic_direction(grid, "grid")
    FT = float(eltype(grid))
    centers = collect(FT, grid)
    edges = Vector{FT}(undef, length(centers) + 1)
    half = FT(0.5)
    @inbounds begin
        edges[1] = centers[1] - half * (centers[2] - centers[1])
        for index in 2:length(centers)
            edges[index] = half * (centers[index - 1] + centers[index])
        end
        edges[end] = centers[end] +
            half * (centers[end] - centers[end - 1])
    end
    return edges
end

function _target_cell_edges(target_grid, configuration)
    inferred = configuration.cell_edges === nothing
    edges = inferred ?
        spectral_cell_edges(target_grid) :
        collect(float(eltype(target_grid)), configuration.cell_edges)
    length(edges) == length(target_grid) + 1 || throw(DimensionMismatch(
        "cell_edges must contain one more element than the target grid",
    ))
    direction = length(target_grid) == 1 ?
        _monotonic_direction(edges, "cell_edges") :
        _monotonic_direction(target_grid, "target grid")
    _monotonic_direction(edges, "cell_edges") == direction ||
        throw(ArgumentError(
            "cell_edges and target grid must have the same ordering",
        ))
    lower_centers = min.(view(edges, 1:length(target_grid)),
                         view(edges, 2:(length(target_grid) + 1)))
    upper_centers = max.(view(edges, 1:length(target_grid)),
                         view(edges, 2:(length(target_grid) + 1)))
    all((lower_centers .<= target_grid) .&
        (target_grid .<= upper_centers)) ||
        throw(ArgumentError("every target-grid point must lie in its cell"))
    return edges
end

function _fine_sampling_grid(target_grid, sampling::AbstractConservativeSampling)
    configuration = sampling_configuration(sampling)
    edges = _target_cell_edges(target_grid, configuration)
    ascending_edges = first(edges) < last(edges) ? edges : reverse(edges)
    FT = float(promote_type(eltype(target_grid), eltype(edges)))
    lower = FT(first(ascending_edges))
    upper = FT(last(ascending_edges))
    target_width = minimum(abs, diff(ascending_edges))
    requested_step = configuration.fine_step === nothing ?
        target_width / FT(configuration.refinement) :
        FT(configuration.fine_step)
    intervals = max(1, ceil(Int, (upper - lower) / requested_step))
    fine_grid = collect(range(lower, upper; length=intervals + 1))
    return fine_grid, edges
end

function _ascending_grid_and_values(source_grid, source_values)
    direction = _monotonic_direction(source_grid, "source grid")
    size(source_values, 1) == length(source_grid) || throw(DimensionMismatch(
        "the first value dimension must match the source grid",
    ))
    grid = collect(float(eltype(source_grid)), source_grid)
    values = Array(source_values)
    return direction > 0 ? (grid, values) :
        (reverse(grid), reverse(values; dims=1))
end

function _cumulative_piecewise_linear(grid, values::AbstractMatrix)
    output_type = float(promote_type(eltype(grid), eltype(values)))
    accumulator_type = widen(output_type)
    cumulative = zeros(accumulator_type, size(values))
    half = accumulator_type(0.5)
    @inbounds for column in axes(values, 2), index in 2:length(grid)
        cumulative[index, column] = cumulative[index - 1, column] +
            half * (values[index - 1, column] + values[index, column]) *
            (grid[index] - grid[index - 1])
    end
    return cumulative
end

@inline function _piecewise_linear_integral_at(
    coordinate,
    grid,
    values,
    cumulative,
    column,
)
    coordinate <= first(grid) && return zero(eltype(cumulative))
    coordinate >= last(grid) && return cumulative[end, column]
    index = searchsortedlast(grid, coordinate)
    fraction = (coordinate - grid[index]) /
        (grid[index + 1] - grid[index])
    endpoint = muladd(
        fraction,
        values[index + 1, column] - values[index, column],
        values[index, column],
    )
    return cumulative[index, column] +
        (values[index, column] + endpoint) *
        (coordinate - grid[index]) / 2
end

"""
    conservative_resample(source_grid, source_values, target_grid;
                          cell_edges=nothing)

Integrate a resolved point-sampled spectrum as piecewise-linear segments and
return its target-cell averages. `source_values` may be a vector or a matrix
whose first dimension is spectral. Values outside the source-grid coverage are
zero. Source and target grids may be ascending or descending.
"""
function conservative_resample(
    source_grid::AbstractVector,
    source_values::AbstractArray,
    target_grid::AbstractVector;
    cell_edges=nothing,
)
    ndims(source_values) in (1, 2) || throw(DimensionMismatch(
        "source_values must be a vector or matrix",
    ))
    size(source_values, 1) == length(source_grid) || throw(DimensionMismatch(
        "the first value dimension must match the source grid",
    ))
    output_type = float(promote_type(
        eltype(source_grid),
        eltype(source_values),
    ))
    isempty(target_grid) && return similar(
        Array(source_values),
        output_type,
        (0, size(source_values)[2:end]...),
    )
    grid, values_array = _ascending_grid_and_values(
        source_grid,
        source_values,
    )
    was_vector = ndims(values_array) == 1
    values = was_vector ? reshape(values_array, :, 1) : values_array
    configuration = ConservativeGridConfiguration(; cell_edges)
    edges = _target_cell_edges(target_grid, configuration)
    target_direction = first(edges) < last(edges) ? 1 : -1
    ascending_edges = target_direction > 0 ? edges : reverse(edges)
    cumulative = _cumulative_piecewise_linear(grid, values)
    output = Matrix{output_type}(
        undef,
        length(target_grid),
        size(values, 2),
    )
    @inbounds for cell in 1:length(target_grid)
        lower = ascending_edges[cell]
        upper = ascending_edges[cell + 1]
        inverse_width = inv(upper - lower)
        for column in axes(values, 2)
            output[cell, column] = (
                _piecewise_linear_integral_at(
                    upper, grid, values, cumulative, column,
                ) -
                _piecewise_linear_integral_at(
                    lower, grid, values, cumulative, column,
                )
            ) * inverse_width
        end
    end
    target_direction < 0 && (output = reverse(output; dims=1))
    return was_vector ? vec(output) : output
end

function _conservative_cross_section(
    fine_grid,
    fine_cross_section,
    target_grid,
    sampling::ConservativeCrossSectionSampling,
    cell_edges,
)
    return conservative_resample(
        fine_grid,
        fine_cross_section,
        target_grid;
        cell_edges,
    )
end

function _conservative_cross_section(
    fine_grid,
    fine_cross_section,
    target_grid,
    sampling::ConservativeTransmissionSampling,
    cell_edges,
)
    column = sampling.column_density_molecules_cm2
    fine_transmission = exp.(-column .* Array(fine_cross_section))
    mean_transmission = conservative_resample(
        fine_grid,
        fine_transmission,
        target_grid;
        cell_edges,
    )
    lower = floatmin(eltype(mean_transmission))
    upper = one(eltype(mean_transmission))
    return -log.(clamp.(mean_transmission, lower, upper)) ./ column
end

_sampling_architecture_result(::AbstractCrossSectionModel, values) = values

function compute_cross_section(
    model::AbstractCrossSectionModel,
    grid::AbstractVector,
    pressure::Real,
    temperature::Real;
    kwargs...,
)
    return compute_cross_section(
        model,
        grid,
        pressure,
        temperature,
        PointSampling();
        kwargs...,
    )
end

function compute_cross_section(
    model::AbstractCrossSectionModel,
    grid::AbstractVector,
    pressure::Real,
    temperature::Real,
    sampling::AbstractConservativeSampling;
    kwargs...,
)
    isempty(grid) && return compute_cross_section(
        model,
        grid,
        pressure,
        temperature,
        PointSampling();
        kwargs...,
    )
    fine_grid, cell_edges = _fine_sampling_grid(grid, sampling)
    fine_cross_section = compute_cross_section(
        model,
        fine_grid,
        pressure,
        temperature,
        PointSampling();
        kwargs...,
    )
    reduced = _conservative_cross_section(
        fine_grid,
        fine_cross_section,
        grid,
        sampling,
        cell_edges,
    )
    return _sampling_architecture_result(model, reduced)
end
