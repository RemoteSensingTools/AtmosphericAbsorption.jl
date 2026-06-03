#=
The Port contract. A line-list backend (HITRAN, ExoMol, …) is a subtype of
`AbstractLineListPort` that implements the three functions below. Whatever the
backend, callers receive a uniform `LineDatabase` + `AbstractPartitionFunction`, so
the compute core is identical across sources.
=#

"""Supertype for line-list data sources (HITRAN, ExoMol, …)."""
abstract type AbstractLineListPort end

"""
    load_lines(port; mol, iso=-1, ν_min, ν_max, min_strength=0.0, FT=Float64) -> LineDatabase{FT}

Load and filter a line list from `port` into the uniform columnar layout. `iso=-1`
selects all isotopologues; lines are kept where `ν_min ≤ ν0 ≤ ν_max` and
`S ≥ min_strength`.
"""
function load_lines end

"""
    partition_function(port, mol, iso) -> AbstractPartitionFunction

Partition function Q(T) for `(mol, iso)` as provided by this source.
"""
function partition_function end

"""
    source_metadata(port, mol, iso) -> SourceMetadata

Reference state and provenance for `(mol, iso)`.
"""
function source_metadata end
