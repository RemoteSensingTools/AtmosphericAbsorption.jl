"""
    LineLists

The columnar `LineDatabase` (struct-of-arrays) and the `AbstractLineListPort`
contract that every data source (HITRAN, ExoMol, …) implements. Backends live in
`ports/`; this module defines only the shared layout and interface.
"""
module LineLists

using Arrow
using ..Constants: T_REF, P_REF
using ..PartitionFunctions: TIPS2017PF

export AbstractLineListPort, LineDatabase, SourceMetadata,
       load_lines, partition_function, source_metadata,
       HitranPort

include("columnar.jl")
include("interface.jl")
include("ports/hitran/isotope_data.jl")
include("ports/hitran/par_parser.jl")
include("ports/hitran/HitranPort.jl")

end # module
