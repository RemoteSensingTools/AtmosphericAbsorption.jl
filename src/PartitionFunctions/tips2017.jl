#=
TIPS-2017 total internal partition sums (Gamache et al. 2017) for HITRAN
isotopologues, loaded from data/tips2017.arrow (a long (mol, iso, T, Q) table). At
load we group it into per-(mol, iso) temperature/partition series and spline Q(T) on
demand.
=#

const _TIPS_PATH = joinpath(@__DIR__, "..", "..", "data", "tips2017.arrow")
include_dependency(_TIPS_PATH)              # recompile if the table is regenerated

# (mol, iso) => (temperatures, partition sums), built once from the Arrow table.
const _TIPS = let t = Arrow.Table(_TIPS_PATH)
    d = Dict{Tuple{Int32,Int32},Tuple{Vector{Float64},Vector{Float64}}}()
    for k in eachindex(t.mol)
        T, Q = get!(() -> (Float64[], Float64[]), d, (t.mol[k], t.iso[k]))
        push!(T, t.T[k]); push!(Q, t.Q[k])
    end
    d
end

"""
    TIPS2017PF()

HITRAN partition functions from TIPS-2017. Iso-aware: implements
`Q_ratio(pf, mol, iso, T, T_ref)` over the whole table; the (mol, iso)-free `Q`
is intentionally undefined since the partition sum depends on the isotopologue.
"""
struct TIPS2017PF <: AbstractPartitionFunction end

@inline function _tips_series(mol::Integer, iso::Integer)
    s = get(_TIPS, (Int32(mol), Int32(iso)), nothing)
    s === nothing && throw(ArgumentError("TIPS-2017 has no data for (mol=$mol, iso=$iso)"))
    return s
end

"""
    Q_ratio(::TIPS2017PF, mol, iso, T, T_ref) -> Float64

Partition-sum ratio Q(T_ref)/Q(T) for `(mol, iso)`, cubic-splined in temperature.
Throws if either `T` or `T_ref` is outside the tabulated range.
"""
function Q_ratio(::TIPS2017PF, mol::Integer, iso::Integer, T, T_ref)
    TT, TQ = _tips_series(mol, iso)
    Tmin, Tmax = first(TT), last(TT)
    for (label, Tq) in (("T", T), ("T_ref", T_ref))
        Tmin ≤ Tq ≤ Tmax || throw(ArgumentError(
            "TIPS-2017: $label=$Tq K outside [$Tmin, $Tmax] for (mol=$mol, iso=$iso)"))
    end
    spline = CubicSpline(TQ, TT)
    return spline(T_ref) / spline(T)
end
