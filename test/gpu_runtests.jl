#=
GPU validation — kept out of the default suite so CI without a GPU stays light.
Run in an environment that has CUDA.jl (or Metal.jl) alongside the package:

    julia --project=<env-with-CUDA> test/gpu_runtests.jl

Asserts the cross-section kernel on the GPU matches the CPU result to (near) machine
precision in the backend's working precision.
=#
using AtmosphericAbsorption
using AtmosphericAbsorption: SpeedDependentVoigt, HartmannTran
using AtmosphericAbsorption.Architectures: GPU, MetalGPU, CPU
using Test

# Import the backend (if installed) and check functionality in SEPARATE top-level
# statements so the functional() call runs in a newer world than the `using`.
const HAS_CUDA  = Base.find_package("CUDA")  !== nothing
const HAS_METAL = Base.find_package("Metal") !== nothing
HAS_CUDA  && @eval using CUDA
HAS_METAL && @eval using Metal

const arch = if HAS_CUDA && CUDA.functional()
    GPU()
elseif HAS_METAL && Metal.functional()
    MetalGPU()
else
    nothing
end

@testset "GPU kernel matches CPU" begin
    if arch === nothing
        @info "No functional CUDA/Metal backend found; nothing to test."
    else
        @info "Testing on $(arch)"
        FT = arch isa MetalGPU ? Float32 : Float64
        db = LineDatabase(; mol = Int32[2, 2], iso = Int32[1, 1], ν0 = FT[1000, 1000.4],
                          S = FT[1e-21, 2e-21], E_lower = FT[120, 250], g_upper = FT[1, 1],
                          γ_air = FT[0.08, 0.08], γ_self = FT[0, 0], n_air = FT[0.7, 0.7],
                          δ_air = FT[-0.005, -0.005], molar_mass = FT[44, 44],
                          γ2_air = FT[0.01, 0.008], δ2_air = FT[0.003, 0.002],
                          νVC = FT[0.015, 0.012], η = FT[0.3, 0.25],
                          meta = SourceMetadata("synthetic", 296.0, 1013.25))
        pf = TabulatedPF(FT[150, 296, 400], FT[280, 300, 330])
        grid = collect(FT, 998:FT(0.01):1002)
        for profile in (Voigt(), Lorentz(), Doppler(), SpeedDependentVoigt(), HartmannTran())
            cpu = LineByLineModel(db, pf; profile, wing_cutoff = FT(40), architecture = CPU())
            gpu = LineByLineModel(db, pf; profile, wing_cutoff = FT(40), architecture = arch)
            σc = compute_cross_section(cpu, grid, 700.0, 260.0)
            σg = Array(compute_cross_section(gpu, grid, 700.0, 260.0))
            @test eltype(σg) === FT
            @test isapprox(σg, σc; rtol = FT === Float64 ? 1e-10 : 1e-4)
        end
    end
end
