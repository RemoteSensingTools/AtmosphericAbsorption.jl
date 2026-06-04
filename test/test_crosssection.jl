using AtmosphericAbsorption
using AtmosphericAbsorption.LineShapes: voigt, lorentz, doppler, HumlicekWeideman32
using AtmosphericAbsorption.Constants: SQRT_LN2, SQRT_2LN2, C_LIGHT, K_BOLTZMANN, AMU
using Test

# Doppler HWHM for a line at ν0 of molar mass M[g/mol] at temperature T — mirrors prepare.
doppler_hwhm(::Type{FT}, ν0, M, T) where {FT} =
    FT(ν0) * FT(SQRT_2LN2) / FT(C_LIGHT) * sqrt(FT(K_BOLTZMANN) * FT(T) / FT(AMU)) / sqrt(FT(M))

# Build a flat-Q partition function (Q_ratio ≡ 1) so T=T_ref gives unit corrections.
flatpf(::Type{FT}) where {FT} = TabulatedPF(FT[100, 300, 500], FT[1, 1, 1])

function oneline_db(::Type{FT}; ν0 = 1000, S = 1e-21, γ_air = 0.1, M = 44,
                    E_lower = 0, n_air = 0, δ_air = 0) where {FT}
    LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = FT[ν0], S = FT[S],
                 E_lower = FT[E_lower], g_upper = FT[1], γ_air = FT[γ_air], γ_self = FT[0],
                 n_air = FT[n_air], δ_air = FT[δ_air], molar_mass = FT[M],
                 meta = SourceMetadata("synthetic", 296.0, 1013.25))
end

@testset "cross-section compute core" begin
    @testset "single line == S·voigt over the window ($FT)" for FT in (Float32, Float64)
        model = LineByLineModel(oneline_db(FT), flatpf(FT); profile = Voigt(), wing_cutoff = FT(40))
        grid  = collect(FT, 990:FT(0.01):1010)
        σ = compute_cross_section(model, grid, 1013.25, 296.0)   # p_ref, T_ref ⇒ S unchanged
        γd = doppler_hwhm(FT, 1000, 44, 296)
        y  = FT(SQRT_LN2) * FT(0.1) / γd
        ref = [FT(1e-21) * voigt(HumlicekWeideman32(), ν - FT(1000), γd, y) for ν in grid]
        @test eltype(σ) === FT
        @test isapprox(σ, ref; rtol = FT === Float64 ? 1e-12 : 1e-5)
    end

    @testset "lines add and respect the wing cutoff ($FT)" for FT in (Float32, Float64)
        db = LineDatabase(; mol = Int32[2, 2], iso = Int32[1, 1], ν0 = FT[1000, 1000.5],
                          S = FT[1e-21, 2e-21], E_lower = FT[0, 0], g_upper = FT[1, 1],
                          γ_air = FT[0.1, 0.1], γ_self = FT[0, 0], n_air = FT[0, 0],
                          δ_air = FT[0, 0], molar_mass = FT[44, 44],
                          meta = SourceMetadata("synthetic", 296.0, 1013.25))
        model = LineByLineModel(db, flatpf(FT); profile = Lorentz(), wing_cutoff = FT(5))
        grid  = collect(FT, 999:FT(0.01):1001)
        σ = compute_cross_section(model, grid, 1013.25, 296.0)
        γl = FT(0.1)
        ref = [FT(1e-21) * lorentz(ν - FT(1000), γl) + FT(2e-21) * lorentz(ν - FT(1000.5), γl) for ν in grid]
        @test isapprox(σ, ref; rtol = FT === Float64 ? 1e-12 : 1e-5)

        # a line 100 cm⁻¹ away with a 5 cm⁻¹ cutoff must contribute nothing
        far = LineByLineModel(oneline_db(FT; ν0 = 1100), flatpf(FT); wing_cutoff = FT(5))
        @test all(iszero, compute_cross_section(far, grid, 1013.25, 296.0))
    end

    @testset "temperature + pressure correction and shift ($FT)" for FT in (Float32, Float64)
        # Non-zero E_lower/n_air/δ_air exercise the Boltzmann + stimulated-emission +
        # width-scaling + pressure-shift paths that the T=T_ref tests leave at unity.
        E, n, δ = 500, 0.75, -0.01
        db = oneline_db(FT; E_lower = E, n_air = n, δ_air = δ)
        T, p = FT(250), FT(2 * 1013.25)
        Tref, pref = FT(296), FT(1013.25)
        model = LineByLineModel(db, flatpf(FT); profile = Voigt(), wing_cutoff = FT(40))
        νs = FT(1000) + (p / pref) * FT(δ)                 # shifted center
        grid = collect(FT, 999:FT(0.002):1001)
        σ = compute_cross_section(model, grid, p, T)
        c₂ = FT(AtmosphericAbsorption.Constants.C2_RAD)
        Scorr = FT(1e-21) * exp(c₂ * FT(E) * (1 / Tref - 1 / T)) *
                (-expm1(-c₂ * FT(1000) / T)) / (-expm1(-c₂ * FT(1000) / Tref))
        γl = FT(0.1) * (p / pref) * (Tref / T)^FT(n)
        γd = doppler_hwhm(FT, 1000, 44, 250)
        y  = FT(SQRT_LN2) * γl / γd
        ip = argmin(abs.(grid .- νs))
        expected = Scorr * voigt(HumlicekWeideman32(), grid[ip] - νs, γd, y)
        @test isapprox(σ[ip], expected; rtol = FT === Float64 ? 1e-10 : 1e-4)
        @test argmax(σ) == ip                              # peak sits at the shifted center
    end

    @testset "self-broadening n_self + δ_self with vmr ($FT)" for FT in (Float32, Float64)
        # A distinct self exponent/shift plus vmr>0 exercise the VMR-weighted width/shift the
        # air-only tests leave at unity. Matches HAPI:
        #   Γ0 = p̃·(γa·(1-v)·(Tr/T)^na + γs·v·(Tr/T)^ns),  Δ0 = p̃·(δa·(1-v) + δs·v).
        γa, γs, na, ns, δa, δs, v = 0.09, 0.12, 0.7, 0.5, -0.008, 0.02, 0.3
        db = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = FT[1000], S = FT[1e-21],
                          E_lower = FT[0], g_upper = FT[1], γ_air = FT[γa], γ_self = FT[γs],
                          n_air = FT[na], δ_air = FT[δa], n_self = FT[ns], δ_self = FT[δs],
                          molar_mass = FT[44], meta = SourceMetadata("synthetic", 296.0, 1013.25))
        T, p, Tref, pref = FT(250), FT(2 * 1013.25), FT(296), FT(1013.25)
        pr = p / pref
        model = LineByLineModel(db, flatpf(FT); profile = Voigt(), wing_cutoff = FT(40), vmr = FT(v))
        Γ0 = pr * (FT(γa) * (1 - FT(v)) * (Tref / T)^FT(na) + FT(γs) * FT(v) * (Tref / T)^FT(ns))
        Δ0 = pr * (FT(δa) * (1 - FT(v)) + FT(δs) * FT(v))
        νs = FT(1000) + Δ0
        grid = collect(FT, 999:FT(0.002):1001)
        σ = compute_cross_section(model, grid, p, T)
        c₂ = FT(AtmosphericAbsorption.Constants.C2_RAD)
        Scorr = FT(1e-21) * (-expm1(-c₂ * FT(1000) / T)) / (-expm1(-c₂ * FT(1000) / Tref))  # E=0 ⇒ Boltz=1
        γd = doppler_hwhm(FT, 1000, 44, 250)
        y  = FT(SQRT_LN2) * Γ0 / γd
        ip = argmin(abs.(grid .- νs))
        @test isapprox(σ[ip], Scorr * voigt(HumlicekWeideman32(), grid[ip] - νs, γd, y);
                       rtol = FT === Float64 ? 1e-10 : 1e-4)
        @test argmax(σ) == ip
        # The self term genuinely matters: it differs from the air-only width here.
        @test !isapprox(Γ0, pr * FT(γa) * (Tref / T)^FT(na); rtol = 1e-3)
    end

    @testset "empty grid returns empty ($FT)" for FT in (Float32, Float64)
        model = LineByLineModel(oneline_db(FT), flatpf(FT))
        @test isempty(compute_cross_section(model, FT[], 1013.25, 296.0))
    end

    @testset "type stability of compute_cross_section ($FT)" for FT in (Float32, Float64)
        model = LineByLineModel(oneline_db(FT), flatpf(FT))
        grid  = collect(FT, 999:FT(0.05):1001)
        @test @inferred(compute_cross_section(model, grid, 1013.25, 296.0)) isa Vector{FT}
    end
end
