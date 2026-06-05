using AtmosphericAbsorption
using Test
using LazyArtifacts

# HAPI golden references and other reference data live in a lazy Pkg artifact (kept out of
# the package tree). Resolve it once; if it can't be fetched (offline), skip the tests that
# need it rather than fail. `GOLDEN`/`COPAR` are shared by the golden test files below.
const GOLDEN = try
    joinpath(artifact"reference_data", "golden")
catch err
    @warn "reference_data artifact unavailable — skipping golden/reference tests" exception = (err, catch_backtrace())
    nothing
end
const COPAR = GOLDEN === nothing ? "" : joinpath(GOLDEN, "co_2100_2200.par")

@testset "AtmosphericAbsorption" begin
    include("test_cpf.jl")
    include("test_profiles.jl")
    include("test_crosssection.jl")
    include("test_advanced_profiles.jl")
    include("test_xsc.jl")
    include("test_interpolation.jl")
    include("test_continuum.jl")
    if GOLDEN !== nothing
        include("test_pcqsdhc.jl")
        include("test_partition.jl")
        include("test_species.jl")
        include("test_hitran_golden.jl")
        include("test_exomol.jl")
    end
end
