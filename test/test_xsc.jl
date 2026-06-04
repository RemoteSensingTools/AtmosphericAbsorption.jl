using AtmosphericAbsorption
using Printf
using Test

# Write a minimal HITRAN-format .xsc with two (T,p) panels (same molecule and region).
function write_xsc(path)
    open(path, "w") do io
        for (T, base) in ((250.0, 1.0), (300.0, 3.0))
            @printf io "%-20s%10.4f%10.4f%7d%7.1f%6.1f\n" "TESTMOL" 800.0 900.0 12 T 760.0
            for j in 0:11
                @printf io "%10.3E" (base * (1 + 0.1j)) * 1e-18
                (j == 11 || (j + 1) % 10 == 0) && print(io, "\n")
            end
        end
    end
end

@testset "tabulated cross-sections (.xsc)" begin
    path = joinpath(mktempdir(), "test.xsc")
    write_xsc(path)

    bands = read_xsc(path)
    @test length(bands) == 2
    @test bands[1].molecule == "TESTMOL"
    @test (bands[1].νmin, bands[1].νmax) == (800.0, 900.0)
    @test length(bands[1].σ) == 12
    @test getfield.(bands, :T) == [250.0, 300.0]
    @test bands[1].p ≈ 760.0 * (101325 / 760 / 100)            # Torr → hPa (≈ 1013.25)
    @test bands[1].σ[1] ≈ 1e-18 && bands[1].σ[end] ≈ 1e-18 * (1 + 0.1 * 11)

    model = load_xsc(path)
    @test model isa TabulatedCrossSection
    p = bands[1].p
    g = collect(800.0:10.0:900.0)
    s250 = compute_cross_section(model, g, p, 250.0)
    s300 = compute_cross_section(model, g, p, 300.0)
    @test s250[1] ≈ 1e-18                                       # exact tabulated node
    @test maximum(abs.(compute_cross_section(model, g, p, 275.0) .- (s250 .+ s300) ./ 2)) < 1e-30
    @test compute_cross_section(model, g, p, 200.0) ≈ s250      # clamp below the T range
    @test compute_cross_section(model, g, p, 400.0) ≈ s300      # clamp above
    @test all(iszero, compute_cross_section(model, [700.0, 950.0], p, 250.0))  # outside coverage

    m32 = load_xsc(path; FT = Float32)
    @test eltype(compute_cross_section(m32, Float32.(g), 1013.25f0, 250.0f0)) === Float32
end
