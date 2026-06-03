#=
Speed benchmark for AtmosphericAbsorption.jl, vs the HAPI reference timed by
benchmark/hapi_reference/time_hapi.py on the same case. Run:

    julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
    julia --project=benchmark benchmark/benchmarks.jl

Times the end-to-end cross-section (pre-pass + kernel) on a real HITRAN CO2 window on
CPU, and on the GPU when CUDA is available. Reports throughput in line·point/s.
=#
using AtmosphericAbsorption
using AtmosphericAbsorption.Architectures: GPU, CPU
using BenchmarkTools

const GOLDEN = joinpath(@__DIR__, "..", "test", "golden")
const CASE   = "co2_6300_6400_p500_T250"   # 4971 lines, 10001 grid points

const HAS_CUDA = Base.find_package("CUDA") !== nothing
HAS_CUDA && @eval using CUDA

function bench(FT)
    db   = load_lines(HitranPort(joinpath(GOLDEN, CASE * ".par")); mol = 2, iso = 1, FT = FT)
    grid = collect(FT, 6300:FT(0.01):6400)
    work = length(db) * length(grid)

    cpu = LineByLineModel(db, TIPS2017PF(); profile = Voigt(), wing_cutoff = FT(40), architecture = CPU())
    tcpu = @belapsed compute_cross_section($cpu, $grid, 500.0, 250.0)
    println("  CPU  $FT: $(round(1e3 * tcpu, digits = 2)) ms  " *
            "($(round(work / tcpu / 1e9, digits = 2)) Gline·pt/s)")

    if HAS_CUDA && CUDA.functional()
        gpu = LineByLineModel(db, TIPS2017PF(); profile = Voigt(), wing_cutoff = FT(40), architecture = GPU())
        tgpu = @belapsed CUDA.@sync compute_cross_section($gpu, $grid, 500.0, 250.0)
        println("  GPU  $FT: $(round(1e3 * tgpu, digits = 2)) ms  " *
                "($(round(work / tgpu / 1e9, digits = 2)) Gline·pt/s, $(round(tcpu / tgpu, digits = 1))× CPU)")
    end
end

println("Case $CASE — $(length(load_lines(HitranPort(joinpath(GOLDEN, CASE * ".par")); mol=2, iso=1))) lines × 10001 points")
for FT in (Float64, Float32)
    bench(FT)
end
