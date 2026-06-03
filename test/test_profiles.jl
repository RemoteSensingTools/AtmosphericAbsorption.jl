using AtmosphericAbsorption: Doppler, Lorentz, Voigt, HumlicekWeideman32, ErfcxCPF,
                             doppler, lorentz, voigt
using AtmosphericAbsorption.LineShapes: evaluate   # internal kernel dispatch helper
using JET

# Trapezoidal integral of a profile over a symmetric grid — used for ∫ϕ dν ≈ 1.
# Domain spans ±60 (3000× the widths) so the slow Lorentzian wings are captured.
trapz(f, νs) = sum((f(a) + f(b)) * (b - a) / 2 for (a, b) in zip(νs, @view νs[2:end]))

@testset "line profiles" begin
    cpf = HumlicekWeideman32()
    ref = ErfcxCPF()                       # accurate reference for the analytic limits

    @testset "normalization ∫ϕ dν ≈ 1 ($FT)" for FT in (Float32, Float64)
        rtol = FT === Float64 ? 1.5e-3 : 3e-3
        γd, γl = FT(0.02), FT(0.02)
        y  = sqrt(log(FT(2))) * γl / γd
        νs = FT.(-60:0.001:60)
        @test isapprox(trapz(Δν -> doppler(Δν, γd), νs), 1; rtol)
        @test isapprox(trapz(Δν -> lorentz(Δν, γl), νs), 1; rtol)
        @test isapprox(trapz(Δν -> voigt(cpf, Δν, γd, y), νs), 1; rtol)
    end

    @testset "Voigt → Doppler as y → 0 ($FT)" for FT in (Float32, Float64)
        γd = FT(0.04)
        y  = FT(1e-6)                       # near-zero Lorentz width
        for Δν in FT.(-0.1:0.01:0.1)        # within ~2.5 Doppler HWHM (not the noise floor)
            @test isapprox(voigt(ref, Δν, γd, y), doppler(Δν, γd); rtol = FT(1e-4))
        end
    end

    @testset "Voigt(0) → Lorentz(0) as y → ∞ ($FT)" for FT in (Float32, Float64)
        # At line center, V(0) → 1/(π γl) = lorentz(0, γl) when Doppler is negligible.
        γl = FT(0.05)
        γd = FT(1e-4)
        y  = sqrt(log(FT(2))) * γl / γd
        @test isapprox(voigt(ref, zero(FT), γd, y), lorentz(zero(FT), γl); rtol = FT(1e-2))
    end

    @testset "evaluate dispatches to the right profile ($FT)" for FT in (Float32, Float64)
        ν0, γd, Γ0, Δ0 = FT(1000), FT(0.04), FT(0.03), FT(-0.01)
        νI = FT(1000.02)
        p  = (γd = γd, Γ0 = Γ0, Γ2 = zero(FT), Δ0 = Δ0,
              Δ2 = zero(FT), νVC = zero(FT), η = zero(FT), Y = zero(FT))
        y  = sqrt(log(FT(2))) * Γ0 / γd
        @test evaluate(Doppler(), cpf, νI, ν0, p) == doppler(νI - (ν0 + Δ0), γd)
        @test evaluate(Lorentz(), cpf, νI, ν0, p) == lorentz(νI - (ν0 + Δ0), Γ0)
        @test evaluate(Voigt(),   cpf, νI, ν0, p) ≈ voigt(cpf, νI - (ν0 + Δ0), γd, y)
        @test @inferred(evaluate(Voigt(), cpf, νI, ν0, p)) isa FT
    end

    let p = (γd = 0.04, Γ0 = 0.03, Γ2 = 0.0, Δ0 = -0.01, Δ2 = 0.0, νVC = 0.0, η = 0.0, Y = 0.0)
        @test_opt evaluate(Voigt(), cpf, 1000.02, 1000.0, p)
    end
end
