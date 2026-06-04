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
using Printf: @sprintf, @printf
import Downloads, SHA, Dates, TOML
using ..Constants: T_REF, P_REF, C_LIGHT, C2_RAD
using ..PartitionFunctions: AbstractPartitionFunction, TIPS2017PF, TIPS2021PF, TabulatedPF, pf_name

export AbstractLineListPort, LineDatabase, SourceMetadata,
       load_lines, partition_function, source_metadata,
       HitranPort, ExoMolPort, fetch_hitran, activate_hitran!,
       fetch_hitran_nonvoigt, load_hitran_nonvoigt,
       molecules, molecule_number, molecule_symbol, isotopologue,
       register_molecule!, register_isotopologue!, resolve_molecule, resolve_isotopologue

include("species.jl")
include("interface.jl")
include("columnar.jl")
include("ports/hitran/isotope_data.jl")
include("ports/hitran/par_parser.jl")
include("ports/hitran/HitranPort.jl")
include("ports/hitran/auth.jl")
include("ports/hitran/download.jl")
include("ports/hitran/nonvoigt.jl")
include("ports/exomol/io.jl")
include("ports/exomol/ExoMolPort.jl")

end # module
