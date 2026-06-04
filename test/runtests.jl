using AtmosphericAbsorption
using Test

@testset "AtmosphericAbsorption" begin
    include("test_cpf.jl")
    include("test_profiles.jl")
    include("test_pcqsdhc.jl")
    include("test_partition.jl")
    include("test_crosssection.jl")
    include("test_advanced_profiles.jl")
    include("test_xsc.jl")
    include("test_species.jl")
    include("test_hitran_golden.jl")
    include("test_exomol.jl")
    include("test_continuum.jl")
end
