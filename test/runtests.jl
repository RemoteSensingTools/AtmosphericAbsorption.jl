using AtmosphericAbsorption
using Test

@testset "AtmosphericAbsorption" begin
    include("test_cpf.jl")
    include("test_profiles.jl")
    include("test_partition.jl")
end
