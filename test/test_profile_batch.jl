# Batched profile queries vs the scalar path they replace.
#
# The batched kernels claim BIT-IDENTICAL results to per-point scalar queries
# on the CPU backend (same expression tree, same order, same bracket/clamp and
# zero-outside semantics). That claim is the test: `==`, not `isapprox`, on
# every element, for both tabulated model types, including edge cases the
# kernels special-case (query on the table's own nu grid, points outside the
# tabulated (p, T, vmr) ranges, a grid straddling the nu range, single-point
# batches, and a multi-column "many cells" batch).

using Test
using AtmosphericAbsorption
using AtmosphericAbsorption: Architectures
using AtmosphericAbsorption.Crosssections: InterpolationModel, AbscoLUT,
                                           compute_cross_section,
                                           compute_cross_section_profile

@testset "batched profile queries" begin
    FT = Float32
    rng_ν = collect(FT, range(6100.0, 6300.0; length = 401))

    # -- synthetic InterpolationModel (regular p/T grid) ----------------------
    ps = collect(FT, [1.0, 10.0, 100.0, 500.0, 1000.0])
    Ts = collect(FT, 180.0:20.0:320.0)
    σim = FT.(reshape(1 .+ sin.(1:length(rng_ν) * length(ps) * length(Ts)),
                      length(rng_ν), length(ps), length(Ts))) .* FT(1e-25)
    im = InterpolationModel(rng_ν, ps, Ts, σim, FT(0.0), Architectures.CPU())

    # -- synthetic AbscoLUT (per-pressure T axis + broadener axis) -----------
    Tmat = Matrix{FT}(undef, 6, length(ps))
    for (j, p) in enumerate(ps)
        Tmat[:, j] = collect(FT, range(170 + j * 5, 320 + j * 5; length = 6))
    end
    vmrs = collect(FT, [0.0, 0.03, 0.06])
    σab = FT.(reshape(1 .+ cos.(1:length(rng_ν) * length(vmrs) * 6 * length(ps)),
                      length(rng_ν), length(vmrs), 6, length(ps))) .* FT(1e-25)
    lut = AbscoLUT(2, 1, rng_ν, ps, Tmat, vmrs, σab; architecture = Architectures.CPU())

    # query points: interior, on-node, below-range, above-range — the clamp
    # paths all get exercised
    qp = FT[0.5, 1.0, 37.0, 100.0, 750.0, 1500.0]
    qT = FT[150.0, 200.0, 261.0, 275.0, 301.0, 400.0]
    qv = FT[0.0, 0.01, 0.03, 0.045, 0.06, 0.2]

    # query grids: off-node, the table's own nu (fast path), and straddling
    grids = Dict(
        "off-node grid"  => collect(FT, range(6100.05, 6299.9; length = 777)),
        "table nu grid"  => rng_ν,
        "straddling"     => collect(FT, range(6050.0, 6350.0; length = 333)),
    )

    for (name, grid) in grids
        @testset "InterpolationModel, $name" begin
            batched = compute_cross_section_profile(im, grid, qp, qT)
            @test size(batched) == (length(grid), length(qp))
            for k in eachindex(qp)
                scalar = compute_cross_section(im, grid, qp[k], qT[k])
                @test batched[:, k] == scalar
            end
        end
        @testset "AbscoLUT, $name" begin
            batched = compute_cross_section_profile(lut, grid, qp, qT; vmr = qv)
            for k in eachindex(qp)
                scalar = compute_cross_section(lut, grid, qp[k], qT[k];
                                               vmr = qv[k], interp = :linear)
                @test batched[:, k] == scalar
            end
        end
    end

    @testset "no broadener vector -> vmr = 0, matching the scalar default" begin
        grid = grids["off-node grid"]
        batched = compute_cross_section_profile(lut, grid, qp, qT)
        for k in eachindex(qp)
            @test batched[:, k] == compute_cross_section(lut, grid, qp[k], qT[k];
                                                         vmr = 0, interp = :linear)
        end
    end

    @testset "multi-cell batch == concatenated single-cell batches" begin
        grid = grids["off-node grid"]
        # two "cells" of 6 layers each, concatenated on the point axis
        p2 = vcat(qp, reverse(qp)); T2 = vcat(qT, reverse(qT)); v2 = vcat(qv, reverse(qv))
        big = compute_cross_section_profile(lut, grid, p2, T2; vmr = v2)
        a = compute_cross_section_profile(lut, grid, qp, qT; vmr = qv)
        b = compute_cross_section_profile(lut, grid, reverse(qp), reverse(qT); vmr = reverse(qv))
        @test big == hcat(a, b)
    end

    @testset "single point and empty batch" begin
        grid = grids["off-node grid"]
        one = compute_cross_section_profile(im, grid, qp[3:3], qT[3:3])
        @test one[:, 1] == compute_cross_section(im, grid, qp[3], qT[3])
        empty = compute_cross_section_profile(im, grid, FT[], FT[])
        @test size(empty) == (length(grid), 0)
    end

    @testset "argument validation" begin
        grid = grids["off-node grid"]
        @test_throws ArgumentError compute_cross_section_profile(im, grid, qp, qT[1:2])
        @test_throws ArgumentError compute_cross_section_profile(im, grid, qp, qT; vmr = qv)
        @test_throws ArgumentError compute_cross_section_profile(lut, grid, qp, qT; vmr = qv[1:2])
    end

    @testset ":cubic falls back to the scalar path, unchanged" begin
        grid = grids["off-node grid"]
        batched = compute_cross_section_profile(lut, grid, qp, qT; vmr = qv, interp = :cubic)
        for k in eachindex(qp)
            @test batched[:, k] == compute_cross_section(lut, grid, qp[k], qT[k];
                                                         vmr = qv[k], interp = :cubic)
        end
    end
end
println("batched profile tests passed")
