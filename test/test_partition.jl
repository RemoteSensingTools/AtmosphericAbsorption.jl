using AtmosphericAbsorption
using AtmosphericAbsorption: TabulatedPF, Q_ratio
import AtmosphericAbsorption: Q_ratio as Q_ratio_ext
using AtmosphericAbsorption.PartitionFunctions: AbstractPartitionFunction
using Test

# Stub iso-aware backend to confirm the 5-arg Q_ratio dispatch reaches an override.
struct StubISOPF <: AbstractPartitionFunction end
Q_ratio_ext(::StubISOPF, mol::Integer, iso::Integer, T, T_ref) = float(iso)

@testset "TabulatedPF" begin
    # Q(T) = T^1.5 ⇒ Q_ratio(T, T_ref) = (T_ref/T)^1.5, recovered through the spline.
    Ts = collect(150.0:1.0:350.0)
    pf = TabulatedPF(Ts, Ts .^ 1.5)
    for (T, Tref) in ((220.0, 296.0), (296.0, 296.0), (310.0, 250.0))
        @test isapprox(Q_ratio(pf, T, Tref), (Tref / T)^1.5; rtol = 1e-6)
    end
    @test Q_ratio(pf, 296.0, 296.0) == 1

    # 5-arg fallback ignores (mol,iso); an iso-aware override is reached when defined.
    @test Q_ratio(pf, 2, 1, 220.0, 296.0) == Q_ratio(pf, 220.0, 296.0)
    @test Q_ratio(StubISOPF(), 2, 7, 250.0, 296.0) == 7.0
end

@testset "TIPS-2021 (latest edition) + 2017 fallback" begin
    p21, p17 = TIPS2021PF(), TIPS2017PF()
    @test Q_ratio(p21, 1, 1, 296.0, 296.0) == 1               # ref/ref = 1
    # Latest edition is the HitranPort default.
    port = HitranPort(joinpath(@__DIR__, "golden", "co_2100_2200.par"))
    @test partition_function(port, 5, 1) isa TIPS2021PF
    # Physical and close to the previous edition for the major isotopologues.
    for (m, i) in ((1, 1), (2, 1), (5, 1), (7, 1))
        r21 = Q_ratio(p21, m, i, 250.0, 296.0)
        @test r21 > 1                                          # Q(296) > Q(250)
        @test isapprox(r21, Q_ratio(p17, m, i, 250.0, 296.0); rtol = 1e-2)
    end
    # TIPS-2021 has no online Q-file for minor O₃ (3,10); it falls back to the 2017 series.
    @test Q_ratio(p21, 3, 10, 250.0, 296.0) == Q_ratio(p17, 3, 10, 250.0, 296.0)
    # Out-of-range temperature throws.
    @test_throws ArgumentError Q_ratio(p21, 1, 1, 6000.0, 296.0)
end

@testset "partition function drives line-strength T-dependence (within a band)" begin
    # The integrated cross-section of one line equals its T-corrected strength S(T), so the
    # area ratio at two temperatures must equal Q(T_ref)/Q(T)·exp(c₂E(1/T_ref-1/T))·stim —
    # i.e. the TIPS partition function actually enters the line strength end-to-end.
    db = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = [2000.0], S = [1e-21],
                      E_lower = [800.0], g_upper = [1.0], γ_air = [0.05], γ_self = [0.0],
                      n_air = [0.7], δ_air = [0.0], molar_mass = [44.0],
                      meta = SourceMetadata("synthetic", 296.0, 1013.25))
    model = LineByLineModel(db, TIPS2017PF(); profile = Voigt(), wing_cutoff = 150.0)
    dν = 0.002
    grid = collect(1850.0:dν:2150.0)
    area(T) = sum(compute_cross_section(model, grid, 1013.25, T)) * dν
    Tref, T = 296.0, 230.0
    Qr = Q_ratio(TIPS2017PF(), 2, 1, T, Tref)                       # Q(296)/Q(230)
    c₂ = AtmosphericAbsorption.Constants.C2_RAD
    expected = Qr * exp(c₂ * 800.0 * (1 / Tref - 1 / T)) *
               (-expm1(-c₂ * 2000.0 / T)) / (-expm1(-c₂ * 2000.0 / Tref))
    @test area(T) / area(Tref) ≈ expected rtol = 5e-3
    @test !(Qr ≈ 1)                                                 # Q genuinely varies with T
end
