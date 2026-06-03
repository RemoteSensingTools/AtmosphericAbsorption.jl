using AtmosphericAbsorption: TabulatedPF, Q_ratio

@testset "TabulatedPF" begin
    # Q(T) = T^1.5 ⇒ Q_ratio(T, T_ref) = (T_ref/T)^1.5, recovered through the spline.
    Ts = collect(150.0:1.0:350.0)
    pf = TabulatedPF(Ts, Ts .^ 1.5)
    for (T, Tref) in ((220.0, 296.0), (296.0, 296.0), (310.0, 250.0))
        @test isapprox(Q_ratio(pf, T, Tref), (Tref / T)^1.5; rtol = 1e-6)
    end
    @test Q_ratio(pf, 296.0, 296.0) == 1
end
