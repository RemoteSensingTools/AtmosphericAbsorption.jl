using AtmosphericAbsorption
using Test

@testset "interpolation model (precomputed LUT)" begin
    db = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = [1000.0], S = [1e-21],
                      E_lower = [300.0], g_upper = [1.0], γ_air = [0.08], γ_self = [0.12],
                      n_air = [0.7], δ_air = [-0.005], molar_mass = [44.0],
                      meta = SourceMetadata("synthetic", 296.0, 1013.25))
    lbl = LineByLineModel(db; profile = Voigt(), wing_cutoff = 40.0)
    ν   = collect(995.0:0.02:1005.0)
    ps  = [200.0, 500.0, 1013.25]
    Ts  = [220.0, 260.0, 300.0]
    im  = build_interpolation_model(lbl, ν, ps, Ts)
    @test im isa InterpolationModel && size(im.σ) == (length(ν), 3, 3)
    @test_throws ArgumentError build_interpolation_model(lbl, ν, reverse(ps), Ts)

    # At a (p, T) node, querying the table reproduces the line-by-line model exactly.
    @test compute_cross_section(im, ν, 500.0, 260.0)    ≈ compute_cross_section(lbl, ν, 500.0, 260.0)
    @test compute_cross_section(im, ν, 1013.25, 300.0)  ≈ compute_cross_section(lbl, ν, 1013.25, 300.0)

    # Between nodes it interpolates: not equal to the exact sum, but the integrated
    # strength (which is p-independent and smooth in T) tracks it closely.
    dν = 0.02
    σi = compute_cross_section(im, ν, 350.0, 240.0)
    σl = compute_cross_section(lbl, ν, 350.0, 240.0)
    @test all(isfinite, σi) && all(≥(0), σi) && σi != σl
    @test isapprox(sum(σi) * dν, sum(σl) * dν; rtol = 0.05)

    # Clamps outside the tabulated p/T range; resamples onto a different grid.
    @test compute_cross_section(im, ν, 50.0, 260.0)  ≈ compute_cross_section(im, ν, 200.0, 260.0)   # p below
    @test compute_cross_section(im, ν, 500.0, 400.0) ≈ compute_cross_section(im, ν, 500.0, 300.0)   # T above
    g2 = collect(996.0:0.5:1004.0)
    @test length(compute_cross_section(im, g2, 500.0, 260.0)) == length(g2)
    # Wavenumbers outside the table's ν range read as zero (no absorption recorded there).
    @test compute_cross_section(im, [900.0, 1100.0], 500.0, 260.0) == zeros(2)

    # vmr is baked at build time: stored on the table, and it threads into the cube
    # (self-broadening makes vmr=0.3 differ from the default vmr=0).
    im_wet = build_interpolation_model(lbl, ν, ps, Ts; vmr = 0.3)
    @test im.vmr == 0.0 && im_wet.vmr == 0.3
    @test compute_cross_section(im_wet, ν, 500.0, 260.0) != compute_cross_section(im, ν, 500.0, 260.0)

    # Single pressure node exercises the degenerate-bracket path.
    im_1p = build_interpolation_model(lbl, ν, [500.0], Ts)
    @test compute_cross_section(im_1p, ν, 500.0, 260.0) ≈ compute_cross_section(lbl, ν, 500.0, 260.0)

    # Persist / restore.
    path = joinpath(mktempdir(), "lut.bin")
    save_interpolation_model(path, im)
    im2 = load_interpolation_model(path)
    @test im2.σ == im.σ && im2.ν == im.ν && im2.p == im.p && im2.T == im.T
    @test compute_cross_section(im2, ν, 500.0, 260.0) == compute_cross_section(im, ν, 500.0, 260.0)

    # Float32 table (the cube's type follows the model's).
    db32 = LineDatabase(; mol = Int32[2], iso = Int32[1], ν0 = Float32[1000], S = Float32[1e-21],
                        E_lower = Float32[300], g_upper = Float32[1], γ_air = Float32[0.08],
                        γ_self = Float32[0], n_air = Float32[0.7], δ_air = Float32[-0.005],
                        molar_mass = Float32[44], meta = SourceMetadata("synthetic", 296.0, 1013.25))
    im32 = build_interpolation_model(LineByLineModel(db32; profile = Voigt(), wing_cutoff = 40.0f0),
                                     Float32.(ν), Float32.(ps), Float32.(Ts))
    @test eltype(im32.σ) === Float32
    @test eltype(compute_cross_section(im32, Float32.(ν), 500.0f0, 260.0f0)) === Float32
end
