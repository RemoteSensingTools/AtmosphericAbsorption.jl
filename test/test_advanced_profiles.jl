using AtmosphericAbsorption
using AtmosphericAbsorption.LineShapes: pcqsdhc, HumlicekWeideman32
using AtmosphericAbsorption.Constants: SQRT_2LN2, C_LIGHT, K_BOLTZMANN, AMU
using Test

γd_of(::Type{FT}, ν0, M, T) where {FT} =
    FT(ν0) * FT(SQRT_2LN2) / FT(C_LIGHT) * sqrt(FT(K_BOLTZMANN) * FT(T) / FT(AMU)) / sqrt(FT(M))

htpf(::Type{FT}) where {FT} = TabulatedPF(FT[150, 296, 400], FT[1, 1, 1])  # flat Q ⇒ Q_ratio≡1

@testset "HT/SDV through the cross-section kernel" begin
    @testset "HartmannTran end-to-end == Σ S·real(pcqsdhc) ($FT)" for FT in (Float32, Float64)
        db = LineDatabase(; mol = Int32[2, 2], iso = Int32[1, 1], ν0 = FT[1000, 1000.3],
                          S = FT[1e-21, 7e-22], E_lower = FT[0, 0], g_upper = FT[1, 1],
                          γ_air = FT[0.03, 0.035], γ_self = FT[0, 0], n_air = FT[0, 0],
                          δ_air = FT[-0.01, -0.008], molar_mass = FT[44, 44],
                          γ2_air = FT[0.006, 0.005], δ2_air = FT[0.002, 0.001],
                          νVC = FT[0.012, 0.01], η = FT[0.3, 0.25],
                          meta = SourceMetadata("synthetic", 296.0, 1013.25))
        cpf = HumlicekWeideman32()
        model = LineByLineModel(db, htpf(FT); profile = HartmannTran(), cpf, wing_cutoff = FT(40))
        grid = collect(FT, 999:FT(0.01):1001)
        σ = compute_cross_section(model, grid, 1013.25, 296.0)   # p_ref, T_ref ⇒ params unscaled
        ref = [sum(db.S[j] * real(pcqsdhc(cpf, db.ν0[j], γd_of(FT, db.ν0[j], 44, 296),
                   db.γ_air[j], db.γ2_air[j], db.δ_air[j], db.δ2_air[j], db.νVC[j], db.η[j], νi))
                   for j in 1:2) for νi in grid]
        @test isapprox(σ, ref; rtol = FT === Float64 ? 1e-10 : 1e-4)
    end

    @testset "line mixing folds in through the kernel ($FT)" for FT in (Float32, Float64)
        db = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = FT[1000], S = FT[1e-21],
                          E_lower = FT[0], g_upper = FT[1], γ_air = FT[0.03], γ_self = FT[0],
                          n_air = FT[0], δ_air = FT[-0.01], molar_mass = FT[44],
                          γ2_air = FT[0.006], δ2_air = FT[0.002], νVC = FT[0.012], η = FT[0.3],
                          Y_LM = FT[0.4], meta = SourceMetadata("synthetic", 296.0, 1013.25))
        cpf = HumlicekWeideman32()
        model = LineByLineModel(db, htpf(FT); profile = HartmannTran(), cpf, wing_cutoff = FT(40))
        grid = collect(FT, 999:FT(0.01):1001)
        σ = compute_cross_section(model, grid, 1013.25, 296.0)   # Y = Y_LM at p_ref, T_ref
        γd = γd_of(FT, 1000, 44, 296)
        ref = [let LS = pcqsdhc(cpf, FT(1000), γd, FT(0.03), FT(0.006), FT(-0.01),
                                FT(0.002), FT(0.012), FT(0.3), νi)
                   FT(1e-21) * (real(LS) + FT(0.4) * imag(LS))
               end for νi in grid]
        @test isapprox(σ, ref; rtol = FT === Float64 ? 1e-10 : 1e-4)
    end

    @testset "Voigt+mixing == HartmannTran+mixing when advanced cols are zero ($FT)" for FT in (Float32, Float64)
        # With only Y_LM set, the fast Voigt path and the pcqsdhc (HT) path must agree.
        db = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = FT[1000], S = FT[1e-21],
                          E_lower = FT[0], g_upper = FT[1], γ_air = FT[0.05], γ_self = FT[0],
                          n_air = FT[0.7], δ_air = FT[-0.01], molar_mass = FT[44],
                          Y_LM = FT[0.35], meta = SourceMetadata("synthetic", 296.0, 1013.25))
        grid = collect(FT, 999:FT(0.01):1001)
        vo = compute_cross_section(LineByLineModel(db, htpf(FT); profile = Voigt(), wing_cutoff = FT(40)), grid, 800.0, 250.0)
        ht = compute_cross_section(LineByLineModel(db, htpf(FT); profile = HartmannTran(), wing_cutoff = FT(40)), grid, 800.0, 250.0)
        @test isapprox(vo, ht; rtol = FT === Float64 ? 1e-6 : 5e-4)
    end

    @testset "HartmannTran reduces to Voigt when advanced cols are zero ($FT)" for FT in (Float32, Float64)
        db = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = FT[1000], S = FT[1e-21],
                          E_lower = FT[0], g_upper = FT[1], γ_air = FT[0.05], γ_self = FT[0],
                          n_air = FT[0.7], δ_air = FT[-0.01], molar_mass = FT[44],
                          meta = SourceMetadata("synthetic", 296.0, 1013.25))
        grid = collect(FT, 999:FT(0.01):1001)
        ht = compute_cross_section(LineByLineModel(db, htpf(FT); profile = HartmannTran(), wing_cutoff = FT(40)), grid, 800.0, 250.0)
        vo = compute_cross_section(LineByLineModel(db, htpf(FT); profile = Voigt(), wing_cutoff = FT(40)), grid, 800.0, 250.0)
        # F32: HT (pcqsdhc) and Voigt call `w` with sign-flipped real args; Re w is even
        # analytically but Weideman32 isn't bit-exactly even, so allow a looser F32 tol.
        @test isapprox(ht, vo; rtol = FT === Float64 ? 1e-6 : 5e-4)
    end
end
