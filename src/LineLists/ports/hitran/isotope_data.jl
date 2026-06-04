#=
HITRAN isotopologue metadata from data/iso_info.arrow: per-(mol, iso) molar mass,
natural abundance, and global id, plus molecule name → id lookup. Loaded once at
module init into a (mol, iso) → row-index map for O(1) access.
=#

const _ISO_PATH  = joinpath(@__DIR__, "..", "..", "..", "..", "data", "iso_info.arrow")
include_dependency(_ISO_PATH)               # recompile if the table is regenerated
const _ISO_TABLE = Arrow.Table(_ISO_PATH)
const _ISO_ROW   = Dict((_ISO_TABLE.mol[k], _ISO_TABLE.iso[k]) => k
                        for k in eachindex(_ISO_TABLE.mol))

@inline function _iso_row(mol::Integer, iso::Integer)
    r = get(_ISO_ROW, (Int32(mol), Int32(iso)), nothing)
    r === nothing && throw(ArgumentError("no isotopologue data for (mol=$mol, iso=$iso)"))
    return r
end

"""
    molar_mass(mol, iso) -> Float64

Molar mass [g/mol] of HITRAN molecule `mol`, isotopologue `iso`.
"""
molar_mass(mol::Integer, iso::Integer) = _ISO_TABLE.molar_mass[_iso_row(mol, iso)]

"""
    abundance(mol, iso) -> Float64

Terrestrial natural abundance of HITRAN molecule `mol`, isotopologue `iso`.
"""
abundance(mol::Integer, iso::Integer) = _ISO_TABLE.abundance[_iso_row(mol, iso)]

"""
    global_id(mol, iso) -> Int

HITRAN global isotopologue id for molecule `mol`, isotopologue `iso`.
"""
global_id(mol::Integer, iso::Integer) = Int(_ISO_TABLE.global_id[_iso_row(mol, iso)])

# molecule_number / molecule_symbol live in species.jl (the backend-agnostic registry).
