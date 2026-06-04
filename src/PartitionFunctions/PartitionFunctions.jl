"""
    PartitionFunctions

Total internal partition sums Q(T) and the temperature-correction ratio
`Q_ratio(pf, T, T_ref) = Q(T_ref)/Q(T)` used to scale line intensities. Backends
implement a uniform interface so HITRAN (TIPS) and ExoMol (.pf) are interchangeable.
`TabulatedPF` covers any (T, Q) table; the TIPS-2017 NetCDF backend lands with the
HITRAN data layer.
"""
module PartitionFunctions

using DataInterpolations: CubicSpline
using Arrow

export AbstractPartitionFunction, TabulatedPF, TIPS2017PF, TIPS2021PF, Q, Q_ratio

"""Supertype for partition-function backends. Implement `Q(pf, T)`; `Q_ratio` follows."""
abstract type AbstractPartitionFunction end

"""
    Q(pf, T)

Total internal partition sum at temperature `T`. The extension point each backend
implements; the default `Q_ratio` is defined in terms of it.
"""
function Q end

"""
    Q_ratio(pf, T, T_ref)

Partition-sum ratio Q(T_ref)/Q(T) — the temperature correction applied to line
intensities. Defined in terms of `Q(pf, T)`; backends may override for speed.
"""
@inline Q_ratio(pf::AbstractPartitionFunction, T, T_ref) = Q(pf, T_ref) / Q(pf, T)

"""
    Q_ratio(pf, mol, iso, T, T_ref)

Per-line ratio used in the pre-pass. The default ignores `(mol, iso)`; iso-aware
backends (e.g. TIPS over a whole molecule) override this method.
"""
@inline Q_ratio(pf::AbstractPartitionFunction, mol::Integer, iso::Integer, T, T_ref) =
    Q_ratio(pf, T, T_ref)

"""
    TabulatedPF(T, Q)

Partition sum from a tabulated `(T, Q)` grid, interpolated with a cubic spline.
Works for any source that ships Q(T) values (e.g. ExoMol `.pf`).
"""
struct TabulatedPF{FT<:AbstractFloat,ITP} <: AbstractPartitionFunction
    T::Vector{FT}
    Qvals::Vector{FT}
    spline::ITP
end

function TabulatedPF(T::AbstractVector, Qvals::AbstractVector)
    FT = float(promote_type(eltype(T), eltype(Qvals)))
    Tf, Qf = collect(FT, T), collect(FT, Qvals)
    return TabulatedPF(Tf, Qf, CubicSpline(Qf, Tf))
end

@inline Q(pf::TabulatedPF, T) = pf.spline(T)

include("tips.jl")

end # module
