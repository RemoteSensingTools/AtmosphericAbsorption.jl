#=
Backend-agnostic species notation. Molecules and isotopologues are selected by readable
symbols (`:CO2`, `:main`) or strings on top of the raw HITRAN integer ids that every port
speaks internally; `:ALL` is the readable spelling of the `-1` "all" sentinel. The maps are
seeded from the canonical HITRAN molecule list + the bundled isotopologue table and are
*overridable* (`register_molecule!`, `register_isotopologue!`), so a project can adopt its own
accepted codes. Ports resolve user selectors through `resolve_molecule`/`resolve_isotopologue`,
so any new port gets the notation for free.
=#

# Canonical HITRAN molecule list (id → formula); seeds the registry, overridable.
const _HITRAN_MOLECULES = [
    (1, "H2O"), (2, "CO2"), (3, "O3"), (4, "N2O"), (5, "CO"), (6, "CH4"), (7, "O2"),
    (8, "NO"), (9, "SO2"), (10, "NO2"), (11, "NH3"), (12, "HNO3"), (13, "OH"), (14, "HF"),
    (15, "HCl"), (16, "HBr"), (17, "HI"), (18, "ClO"), (19, "OCS"), (20, "H2CO"),
    (21, "HOCl"), (22, "N2"), (23, "HCN"), (24, "CH3Cl"), (25, "H2O2"), (26, "C2H2"),
    (27, "C2H6"), (28, "PH3"), (29, "COF2"), (30, "SF6"), (31, "H2S"), (32, "HCOOH"),
    (33, "HO2"), (34, "O"), (35, "ClONO2"), (36, "NOp"), (37, "HOBr"), (38, "C2H4"),
    (39, "CH3OH"), (40, "CH3Br"), (41, "CH3CN"), (42, "CF4"), (43, "C4H2"), (44, "HC3N"),
    (45, "H2"), (46, "CS"), (47, "SO3"), (48, "C2N2"), (49, "COCl2"), (50, "SO"),
    (51, "CH3F"), (52, "GeH4"), (53, "CS2"), (54, "CH3I"), (55, "NF3"), (56, "H3p"),
    (57, "CH3"), (58, "S2"), (59, "COFCl"), (60, "HONO"), (61, "ClNO2"),
]

const MOLECULE_IDS  = Dict{Symbol,Int}(Symbol(n) => i for (i, n) in _HITRAN_MOLECULES)
const _MOLECULE_SYM = Dict{Int,Symbol}(i => Symbol(n) for (i, n) in _HITRAN_MOLECULES)

# Named isotopologue codes (mol, code) => local id, seeded with :main = 1; overridable.
const ISOTOPOLOGUE_CODES = Dict{Tuple{Int,Symbol},Int}((i, :main) => 1 for (i, _) in _HITRAN_MOLECULES)

# Isotopologue chemical formulas for display (e.g. (2, 1) → "(12C)(16O)2"), from the bundled
# HAPI-derived table; user-editable.
const _ISOFORM_PATH = joinpath(@__DIR__, "..", "..", "data", "hitran_isotopologues.csv")
include_dependency(_ISOFORM_PATH)
const _ISOFORM = let d = Dict{Tuple{Int,Int},String}()
    for ln in eachline(_ISOFORM_PATH)
        startswith(ln, "mol") && continue
        a = split(strip(ln), ',')
        length(a) < 3 && continue
        d[(parse(Int, a[1]), parse(Int, a[2]))] = String(strip(a[3]))
    end
    d
end

"""
    register_molecule!(sym, id)

Register or override the molecule symbol `sym` ↦ HITRAN molecule id `id` (e.g.
`register_molecule!(:HDO, 1)`). Affects `resolve_molecule`/`molecule_symbol` thereafter.
"""
register_molecule!(sym::Symbol, id::Integer) =
    (MOLECULE_IDS[sym] = Int(id); _MOLECULE_SYM[Int(id)] = sym; sym)

"""
    register_isotopologue!(mol, code, local_iso)

Register or override a named isotopologue `code` for molecule `mol` (id or symbol) ↦ HITRAN
local isotopologue id `local_iso` (e.g. `register_isotopologue!(:CO2, Symbol("636"), 2)`).
"""
register_isotopologue!(mol, code::Symbol, local_iso::Integer) =
    (ISOTOPOLOGUE_CODES[(resolve_molecule(mol), code)] = Int(local_iso))

"""
    molecule_number(m) -> Int

HITRAN molecule id for a symbol/string formula (`:CO2`/`"CO2"` → 2) or integer (passthrough).
"""
molecule_number(m::Integer) = Int(m)
molecule_number(m::AbstractString) = molecule_number(Symbol(m))
function molecule_number(m::Symbol)
    id = get(MOLECULE_IDS, m, nothing)
    id === nothing && throw(ArgumentError("unknown molecule :$m — register_molecule!(:$m, id) to add it"))
    return id
end

"""
    molecule_symbol(id) -> Symbol

Formula symbol for a HITRAN molecule id (`2` → `:CO2`). Inverse of `molecule_number`.
"""
function molecule_symbol(id::Integer)
    s = get(_MOLECULE_SYM, Int(id), nothing)
    s === nothing && throw(ArgumentError("unknown molecule id $id"))
    return s
end

"""
    isotopologue(mol, iso) -> String

Chemical formula of HITRAN isotopologue `(mol, iso)` (e.g. `(2, 1)` → `"(12C)(16O)2"`); falls
back to `"iso N"` when not tabulated.
"""
isotopologue(mol::Integer, iso::Integer) = get(_ISOFORM, (Int(mol), Int(iso)), "iso $iso")

"""
    resolve_molecule(m) -> Int

Resolve a molecule selector to a HITRAN id: integer passthrough, `:ALL` → -1 ("all"),
symbol/string formula via the registry.
"""
resolve_molecule(m::Integer) = Int(m)
resolve_molecule(m::AbstractString) = molecule_number(m)
resolve_molecule(m::Symbol) = m === :ALL ? -1 : molecule_number(m)

"""
    resolve_isotopologue(mol, i) -> Int

Resolve an isotopologue selector to a HITRAN local id: integer passthrough, `:ALL` → -1,
`:main` → 1, or a registered named code for `mol` (which may itself be a symbol/string/id).
"""
resolve_isotopologue(mol, i::Integer) = Int(i)
function resolve_isotopologue(mol, i::Symbol)
    i === :ALL && return -1
    m  = resolve_molecule(mol)
    id = get(ISOTOPOLOGUE_CODES, (m, i), nothing)
    id === nothing && throw(ArgumentError(
        "unknown isotopologue :$i for molecule $(m) — use :ALL, :main, an id, or register_isotopologue!"))
    return id
end
