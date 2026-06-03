using AtmosphericAbsorption
using Test

const CFIX = joinpath(@__DIR__, "fixtures", "continuum")

@testset "CIA (collision-induced absorption)" begin
    grid = collect(7800.0:1.0:7810.0)
    t = load_cia(joinpath(CFIX, "o2o2_toy.cia"), grid)
    @test t.pair == "O2-O2" && (t.species_a, t.species_b) == ("O2", "O2")
    @test t.block_T == [250.0, 290.0]
    @test cia_cross_section(t, 250.0)[6] == 2.0e-46        # ν=7805, exact block T
    @test cia_cross_section(t, 290.0)[6] == 1.6e-46
    @test cia_cross_section(t, 270.0)[6] ≈ 1.8e-46         # T midpoint
    @test cia_cross_section(t, 200.0)[6] == 2.0e-46        # below range → constant
    @test cia_cross_section(t, 350.0)[6] == 1.6e-46        # above range → constant
    # a window outside every band returns zeros, not an error
    t2 = load_cia(joinpath(CFIX, "o2o2_toy.cia"), collect(100.0:1.0:200.0))
    @test all(iszero, cia_cross_section(t2, 250.0))

    # multi-band window: each band interpolates in T over ITS OWN temperature set.
    # band1 (7800-7810) measured at 250/290; band2 (13000-13010) at 220/280.
    g  = 7800.0:5.0:13010.0
    mb = load_cia(joinpath(CFIX, "o2o2_toy.cia"), collect(g))
    σ  = cia_cross_section(mb, 250.0)
    @test σ[findfirst(≈(7805.0), g)] ≈ 2.0e-46     # band1: exact T=250
    @test σ[findfirst(≈(13005.0), g)] ≈ 4.0e-46    # band2: 250 ∈ (220,280) ⇒ midpoint of 5e-46,3e-46
end

@testset "MT_CKD water-vapor continuum" begin
    tbl = load_mtckd()                      # bundled AER v4.2
    @test tbl.T_ref == 296.0 && length(tbl.ν) > 1000
    grid = collect(2400.0:1.0:2600.0)       # 4 µm window
    band = build_mtckd_band(tbl, grid)
    σ = h2o_continuum(band, 296.0, 20.0, 993.0)            # T[K], p_h2o, p_dry [hPa]
    @test all(isfinite, σ) && all(>=(0), σ) && maximum(σ) > 0
    @test maximum(h2o_continuum(band, 250.0, 20.0, 993.0)) > maximum(σ)   # self grows as T↓
    # foreign continuum is linear in dry pressure
    σ1 = h2o_continuum(band, 296.0, 0.0, 993.0)
    σ2 = h2o_continuum(band, 296.0, 0.0, 1986.0)
    @test maximum(σ2) ≈ 2 * maximum(σ1) rtol = 1e-10
    # outside the table range (UV/Vis) → zero
    far = build_mtckd_band(tbl, collect(25000.0:10.0:26000.0))
    @test all(iszero, h2o_continuum(far, 296.0, 20.0, 993.0))
end
