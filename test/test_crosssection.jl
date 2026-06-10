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

    @testset "per-call vmr overrides the model's ($FT)" for FT in (Float32, Float64)
        db = LineDatabase(; mol = Int32[1], iso = Int32[1], ν0 = FT[1000], S = FT[1e-21],
            E_lower = FT[0], g_upper = FT[1], γ_air = FT[0.09], γ_self = FT[0.15], n_air = FT[0.7],
            δ_air = FT[-0.008], n_self = FT[0.5], δ_self = FT[0.02], molar_mass = FT[18],
            meta = SourceMetadata("synthetic", 296.0, 1013.25))
        T, p = FT(260), FT(800)
        grid = collect(FT, 999:FT(0.004):1001)
        m0 = LineByLineModel(db, flatpf(FT); profile = Voigt(), wing_cutoff = FT(40), vmr = FT(0))
        mv = LineByLineModel(db, flatpf(FT); profile = Voigt(), wing_cutoff = FT(40), vmr = FT(0.4))
        # Overriding vmr at call time equals building the model with that vmr…
        @test compute_cross_section(m0, grid, p, T; vmr = FT(0.4)) ≈ compute_cross_section(mv, grid, p, T)
        # …the default uses the model's own vmr, and self-broadening genuinely moves σ.
        @test compute_cross_section(mv, grid, p, T) ≈ compute_cross_section(mv, grid, p, T; vmr = FT(0.4))
        @test !(compute_cross_section(m0, grid, p, T) ≈ compute_cross_section(m0, grid, p, T; vmr = FT(0.4)))
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

    @testset "subset is type-stable and partition-preserving ($FT)" for FT in (Float32, Float64)
        db = oneline_db(FT)                              # mol = 2 (CO2), default TIPS2021PF
        @test @inferred(db[db.mol .== 2]) isa typeof(db)
        @test molecules(db) == [:CO2]
        @test db[db.mol .== 2].partition === db.partition
    end

    @testset "wavelength grid: equivalence to reverse-mapped wavenumber ($FT)" for FT in (Float32, Float64)
        nm_per_m = FT(AtmosphericAbsorption.Constants.NM_PER_M)
        # Put a CO2 line in the band: ~12980 cm⁻¹ ≈ 770.4 nm.
        ν0 = nm_per_m / FT(770.4)
        model = LineByLineModel(oneline_db(FT; ν0 = ν0), flatpf(FT); profile = Voigt(), wing_cutoff = FT(40))
        λ  = collect(FT, 770:FT(0.002):771)                       # nm, ascending
        ν  = nm_per_m ./ λ                                        # cm⁻¹, descending (nm-ascending)
        σλ = compute_cross_section(model, λ, 1013.25, 296.0; wavelength_flag = true)
        # Equivalent wavenumber path: the wavenumber API needs an ascending grid, so sort the
        # cm⁻¹ values, compute, then reverse-map σ back to the nm (input) order via the same
        # permutation. This must reproduce the wavelength-flag result bit-for-bit.
        perm = sortperm(ν)                                       # ascending-cm⁻¹ ordering
        σsorted = compute_cross_section(model, ν[perm], 1013.25, 296.0)
        σν = similar(σsorted); σν[perm] = σsorted                # reverse-map to nm order
        @test σλ == σν                                           # bit-identical (pure reordering)
        @test eltype(σλ) === FT
    end

    @testset "wavelength grid: wing cutoff is windowed in cm⁻¹, not nm ($FT)" for FT in (Float32, Float64)
        nm_per_m = FT(AtmosphericAbsorption.Constants.NM_PER_M)
        λ  = collect(FT, 770:FT(0.001):771)                       # nm → [12970.2, 12987.0] cm⁻¹
        νmin = nm_per_m / FT(771)                                 # low end of the cm⁻¹ window
        # A line 10 cm⁻¹ below the window with a 5 cm⁻¹ cutoff → exactly zero contribution.
        # Its nm-equivalent (≈771.9 nm) sits inside a naive nm-margin window [765, 776], so a
        # buggy "cm⁻¹ cutoff applied to nm values" implementation would wrongly include it.
        ν_far = νmin - FT(5) - FT(10)
        far = LineByLineModel(oneline_db(FT; ν0 = ν_far, γ_air = FT(0.1)), flatpf(FT);
                              profile = Lorentz(), wing_cutoff = FT(5))
        σ = compute_cross_section(far, λ, 1013.25, 296.0; wavelength_flag = true)
        @test all(iszero, σ)                                     # out-of-window ⇒ contributes nothing

        # Sanity: a line *inside* the window does contribute, so the zero above is meaningful.
        ν_in = nm_per_m / FT(770.5)
        near = LineByLineModel(oneline_db(FT; ν0 = ν_in, γ_air = FT(0.1)), flatpf(FT);
                               profile = Lorentz(), wing_cutoff = FT(5))
        @test any(!iszero, compute_cross_section(near, λ, 1013.25, 296.0; wavelength_flag = true))
    end

    @testset "wavelength grid: returned σ aligns with input nm order ($FT)" for FT in (Float32, Float64)
        nm_per_m = FT(AtmosphericAbsorption.Constants.NM_PER_M)
        ν0 = nm_per_m / FT(770.5)
        model = LineByLineModel(oneline_db(FT; ν0 = ν0), flatpf(FT); profile = Voigt(), wing_cutoff = FT(40))
        λ  = collect(FT, 770:FT(0.002):771)                       # nm, ascending
        σ  = compute_cross_section(model, λ, 1013.25, 296.0; wavelength_flag = true)
        @test length(σ) == length(λ)
        # σ aligns element-for-element with the nm grid: the peak sits at the grid point
        # nearest the line's wavelength (770.5 nm), regardless of internal cm⁻¹ reordering.
        @test λ[argmax(σ)] ≈ FT(770.5) atol = FT(0.01)
        # And the same call on a deliberately shuffled nm grid scatters σ back to match it.
        perm = [3, 1, length(λ), 2]                              # an arbitrary reordering of a few points
        λp = λ[perm]
        σp = compute_cross_section(model, λp, 1013.25, 296.0; wavelength_flag = true)
        @test σp == σ[perm]
    end
end
