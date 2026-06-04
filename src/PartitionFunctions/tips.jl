#=
HITRAN total internal partition sums (TIPS), loaded from long (mol, iso, T, Q) Arrow
tables and splined per (mol, iso) on demand. Two editions ship:
  • TIPS2021PF — the latest (Gamache et al. 2021), the HITRAN default, from
    data/tips2021.arrow (built from hitran.org/data/Q by data/_build_tips2021.jl).
  • TIPS2017PF — the previous edition (Gamache et al. 2017), retained because the HAPI
    reference uses it (so cross-validation stays apples-to-apples) and for the handful
    of (mostly minor O₃) isotopologues not yet in HITRAN's online Q-file set.
TIPS2021PF falls back to the 2017 series for any (mol, iso) absent from the 2021 table.
=#

const _TIPS2017_PATH = joinpath(@__DIR__, "..", "..", "data", "tips2017.arrow")
const _TIPS2021_PATH = joinpath(@__DIR__, "..", "..", "data", "tips2021.arrow")
include_dependency(_TIPS2017_PATH)               # recompile if a table is regenerated
include_dependency(_TIPS2021_PATH)

# Group a long (mol, iso, T, Q) Arrow table into (mol, iso) => (temperatures, sums).
function _load_tips(path)
    t = Arrow.Table(path)
    d = Dict{Tuple{Int32,Int32},Tuple{Vector{Float64},Vector{Float64}}}()
    for k in eachindex(t.mol)
        T, Q = get!(() -> (Float64[], Float64[]), d, (t.mol[k], t.iso[k]))
        push!(T, t.T[k]); push!(Q, t.Q[k])
    end
    d
end

const _TIPS2017 = _load_tips(_TIPS2017_PATH)
const _TIPS2021 = _load_tips(_TIPS2021_PATH)

"""
    TIPS2017PF()

HITRAN partition functions from TIPS-2017 (Gamache et al. 2017). Iso-aware via
`Q_ratio(pf, mol, iso, T, T_ref)`. Retained mainly for HAPI cross-validation; new work
should prefer the latest edition, [`TIPS2021PF`](@ref).
"""
struct TIPS2017PF <: AbstractPartitionFunction end

"""
    TIPS2021PF()

HITRAN partition functions from the latest edition, TIPS-2021 (Gamache et al. 2021) —
the `HitranPort` default. Iso-aware via `Q_ratio(pf, mol, iso, T, T_ref)`. For the few
isotopologues not yet in HITRAN's online Q-files it falls back to TIPS-2017.
"""
struct TIPS2021PF <: AbstractPartitionFunction end

const _TIPSPF = Union{TIPS2017PF,TIPS2021PF}

pf_name(::TIPS2017PF) = "TIPS-2017"
pf_name(::TIPS2021PF) = "TIPS-2021"

function _tips_series(::TIPS2017PF, mol::Integer, iso::Integer)
    s = get(_TIPS2017, (Int32(mol), Int32(iso)), nothing)
    s === nothing && throw(ArgumentError("TIPS-2017 has no data for (mol=$mol, iso=$iso)"))
    return s
end

# TIPS-2021 falls back to the 2017 series for isotopologues without an online Q-file.
function _tips_series(::TIPS2021PF, mol::Integer, iso::Integer)
    key = (Int32(mol), Int32(iso))
    s = get(_TIPS2021, key, nothing)
    s === nothing || return s
    s17 = get(_TIPS2017, key, nothing)
    s17 === nothing && throw(ArgumentError("TIPS-2021 has no data for (mol=$mol, iso=$iso)"))
    return s17
end

"""
    Q_ratio(pf::Union{TIPS2017PF,TIPS2021PF}, mol, iso, T, T_ref) -> Float64

Partition-sum ratio Q(T_ref)/Q(T) for `(mol, iso)`, cubic-splined in temperature.
Throws if either `T` or `T_ref` is outside the tabulated range.
"""
function Q_ratio(pf::_TIPSPF, mol::Integer, iso::Integer, T, T_ref)
    TT, TQ = _tips_series(pf, mol, iso)
    Tmin, Tmax = first(TT), last(TT)
    for (label, Tq) in (("T", T), ("T_ref", T_ref))
        Tmin ≤ Tq ≤ Tmax || throw(ArgumentError(
            "$(pf_name(pf)): $label=$Tq K outside [$Tmin, $Tmax] for (mol=$mol, iso=$iso)"))
    end
    spline = CubicSpline(TQ, TT)
    return spline(T_ref) / spline(T)
end
