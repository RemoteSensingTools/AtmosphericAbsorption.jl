#=
The Port contract. A line-list backend (HITRAN, ExoMol, …) is a subtype of
`AbstractLineListPort` that implements the three functions below — the whole story for
adding a new source, much like adding a distribution to Distributions.jl. Whatever the
backend, callers receive a uniform `LineDatabase` (with its partition function attached)
and select species with the same generic notation, so the compute core and the user-facing
API are identical across sources.

To add a port:
  1. `struct MyPort <: AbstractLineListPort … end`.
  2. `load_lines(::MyPort; mol, iso, …)` — resolve the user's selectors through the shared
     `resolve_molecule`/`resolve_isotopologue` (so `:CO2`/`"CO2"`/2/`:ALL` all work), and
     return a `LineDatabase` built with `partition = partition_function(port, mol, iso)`.
  3. `partition_function(::MyPort, mol, iso)` and `source_metadata(::MyPort, mol, iso)`.
The output `LineDatabase` is keyed by HITRAN `(mol, iso)` integer ids (the lingua franca);
a backend with its own naming maps onto those ids (see `ExoMolPort`'s `hitran_mol/iso`).
=#

"""Supertype for line-list data sources (HITRAN, ExoMol, …). See the Port contract above."""
abstract type AbstractLineListPort end

"""
    load_lines(port; mol=:ALL, iso=:ALL, ν_min, ν_max, min_strength=0.0, FT=Float64) -> LineDatabase{FT}

Load and filter a line list from `port` into the uniform columnar layout, **with the
source's partition function attached** (`db.partition`). `mol`/`iso` accept the generic
notation (a symbol `:CO2`/`:main`, a string, an integer id, or `:ALL` for "all"); lines are
kept where `ν_min ≤ ν0 ≤ ν_max` and `S ≥ min_strength`. Implementations resolve selectors via
`resolve_molecule`/`resolve_isotopologue`.
"""
function load_lines end

"""
    partition_function(port, mol, iso) -> AbstractPartitionFunction

Partition function Q(T) for `(mol, iso)` as provided by this source. `load_lines` attaches
its result to the returned `LineDatabase`, so callers rarely need it directly.
"""
function partition_function end

"""
    source_metadata(port, mol, iso) -> SourceMetadata

Reference state and provenance for `(mol, iso)`.
"""
function source_metadata end
