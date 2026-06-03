"""
    LineLists

The columnar `LineDatabase` (struct-of-arrays) and the `AbstractLineListPort`
contract that every data source (HITRAN, ExoMol, …) implements. Backends live in
`ports/`; this module defines only the shared layout and interface.
"""
module LineLists

using Arrow
using CodecBzip2: Bzip2DecompressorStream
using DataInterpolations: CubicSpline
using Scratch: @get_scratch!
using Printf: @sprintf
import Downloads, SHA, Dates, TOML
using ..Constants: T_REF, P_REF, C_LIGHT, C2_RAD
using ..PartitionFunctions: TIPS2017PF, TabulatedPF

export AbstractLineListPort, LineDatabase, SourceMetadata,
       load_lines, partition_function, source_metadata,
       HitranPort, ExoMolPort, fetch_hitran

include("columnar.jl")
include("interface.jl")
include("ports/hitran/isotope_data.jl")
include("ports/hitran/par_parser.jl")
include("ports/hitran/HitranPort.jl")
include("ports/hitran/download.jl")
include("ports/exomol/io.jl")
include("ports/exomol/ExoMolPort.jl")

end # module
