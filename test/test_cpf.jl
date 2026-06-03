using AtmosphericAbsorption: HumlicekWeideman32, ErfcxCPF, w
using JET

# `ErfcxCPF` is the gold reference; `HumlicekWeideman32` is the production strategy.
@testset "complex probability function w(z)" begin
    ref  = ErfcxCPF()
    fast = HumlicekWeideman32()

    @testset "HumlicekWeideman32 ≈ erfcx reference ($FT)" for FT in (Float32, Float64)
        rtol = FT === Float64 ? 1e-4 : 1e-3
        xs = FT.(-30:0.37:30)
        ys = FT.(10.0 .^ (-3:0.5:1.5))      # spans Weideman (small) and Humlíček (large) regions
        @test all(
            isapprox(w(fast, x + im * y), w(ref, x + im * y); rtol)
            for x in xs, y in ys
        )
    end

    @testset "w is type stable and FT-preserving ($FT)" for FT in (Float32, Float64)
        z = FT(1.2) + im * FT(0.3)
        @test @inferred(w(fast, z)) isa Complex{FT}
        @test @inferred(w(ref, z))  isa Complex{FT}
    end

    @test_opt w(fast, 1.2 + 0.3im)
end
