"""
    LineLists

The columnar `LineDatabase` (struct-of-arrays) and the `AbstractLineListPort`
contract that every data source (HITRAN, ExoMol, …) implements. Backends live in
`ports/`; this module defines only the shared layout and interface.
"""
module LineLists

export AbstractLineListPort, LineDatabase, SourceMetadata,
       load_lines, partition_function, source_metadata

include("columnar.jl")
include("interface.jl")

end # module
