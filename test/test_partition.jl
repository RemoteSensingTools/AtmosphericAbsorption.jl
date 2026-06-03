using AtmosphericAbsorption: TabulatedPF, Q_ratio
import AtmosphericAbsorption: Q_ratio as Q_ratio_ext
using AtmosphericAbsorption.PartitionFunctions: AbstractPartitionFunction

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
