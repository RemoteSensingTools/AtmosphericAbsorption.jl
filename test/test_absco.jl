using AtmosphericAbsorption
using Test
using NCDatasets   # enables read_absco (the ABSCO .hdf reader extension)

# Build an AbscoLUT from a known analytic σ on an ABSCO-shaped native grid (per-pressure T axis that
# slides with p, plus an H₂O broadener axis), then check we recover the original spectra exactly at the
# original nodes — the core "convert but stay faithful" requirement — and compare linear vs cubic-in-T.
@testset "AbscoLUT — native-grid recovery + interpolation" begin
    ν   = collect(4000.0:1.0:4050.0)                       # 51 wavenumber nodes
    p   = [100.0, 500.0, 1000.0]                           # ascending pressures [hPa]
    nT  = 5
    Tmn = [180.0, 190.0, 200.0]                            # T_min slides with pressure (ABSCO style)
    T   = Float64[Tmn[ip] + 10 * (it - 1) for it in 1:nT, ip in 1:3]   # T[iT, ip], 10 K apart per p
    vmr = [0.0, 0.03, 0.06]                                # H₂O broadener nodes

    # Separable truth: linear in vmr and p (so those blends are exact), quadratic in T, Gaussian in ν.
    σtrue(x, v, t, pp) = exp(-((x - 4025) / 8)^2) * (1 + 2v) *
                         (1 + 5e-3 * (t - 200) + 2e-5 * (t - 200)^2) * (1 + 5e-4 * (pp - 500))
    σ = Float64[σtrue(ν[i], vmr[kv], T[it, ip], p[ip])
                for i in eachindex(ν), kv in 1:3, it in 1:nT, ip in 1:3]

    lut = AbscoLUT(2, -1, ν, p, T, vmr, σ)                 # CPU (no CUDA loaded)
    @test lut isa AbscoLUT && size(lut.σ) == (length(ν), 3, nT, 3)

    # 1. Exact recovery at every original (p, T, vmr) node — for BOTH interpolation modes.
    for ip in 1:3, it in 1:nT, kv in 1:3, interp in (:linear, :cubic)
        σq = compute_cross_section(lut, ν, p[ip], T[it, ip]; vmr = vmr[kv], interp)
        @test σq ≈ σ[:, kv, it, ip]
    end

    # 2. Off-node: cubic reproduces the quadratic T-dependence (Catmull-Rom is exact for quadratics on
    #    interior intervals), so it beats linear; linear is still reasonable.
    pq, Tq, vq = 500.0, 205.0, 0.03                        # p,vmr on nodes; T interior between 200 & 210
    truth = [σtrue(x, vq, Tq, pq) for x in ν]
    σlin  = compute_cross_section(lut, ν, pq, Tq; vmr = vq, interp = :linear)
    σcub  = compute_cross_section(lut, ν, pq, Tq; vmr = vq, interp = :cubic)
    @test maximum(abs.(σcub .- truth)) < maximum(abs.(σlin .- truth))
    @test isapprox(σlin, truth; rtol = 0.05)

    # 3. Clamps p/T/vmr to the table; ν outside the band reads zero; off-grid resample keeps length.
    @test compute_cross_section(lut, ν, 1.0, 200.0)   ≈ compute_cross_section(lut, ν, p[1], 200.0)   # p below
    @test compute_cross_section(lut, ν, 500.0, 999.0) ≈ compute_cross_section(lut, ν, 500.0, T[end, 2])  # T above
    @test compute_cross_section(lut, [3000.0, 5000.0], 500.0, 200.0) == zeros(2)
    g2 = collect(4010.0:2.0:4040.0)
    @test length(compute_cross_section(lut, g2, 500.0, 200.0)) == length(g2)

    # 4. The H₂O broadener really is a third axis: dry ≠ wet, and the default is dry (vmr=0).
    @test compute_cross_section(lut, ν, 500.0, 200.0; vmr = 0.0) !=
          compute_cross_section(lut, ν, 500.0, 200.0; vmr = 0.06)
    @test compute_cross_section(lut, ν, 500.0, 200.0) ≈ compute_cross_section(lut, ν, 500.0, 200.0; vmr = 0.0)

    # 5. Persist / restore.
    path = joinpath(mktempdir(), "absco.bin")
    save_absco_lut(path, lut)
    lut2 = load_absco_lut(path)
    @test lut2.σ == lut.σ && lut2.T == lut.T && lut2.vmr == lut.vmr
    @test compute_cross_section(lut2, ν, 500.0, 205.0; vmr = 0.03) ==
          compute_cross_section(lut, ν, 500.0, 205.0; vmr = 0.03)

    # 6. Float32 table stays Float32 end-to-end.
    lut32 = AbscoLUT(2, -1, Float32.(ν), Float32.(p), Float32.(T), Float32.(vmr), Float32.(σ))
    @test eltype(lut32.σ) === Float32
    @test eltype(compute_cross_section(lut32, Float32.(ν), 500.0f0, 205.0f0; vmr = 0.03f0)) === Float32
end

# Round-trip through the NCDatasets reader: write a tiny ABSCO-shaped HDF5/NetCDF file with the real
# variable names + dimension order, read it back, and confirm faithful recovery at the original nodes.
@testset "read_absco (synthetic ABSCO fixture)" begin
    nν, nb, nT, nP = 8, 3, 4, 3
    νf = collect(4000.0:1.0:4007.0)
    pf = [100.0, 500.0, 1000.0]                                   # hPa (written as Pa below)
    vf = [0.0, 0.03, 0.06]
    Tf = Float64[170 + 10 * (it - 1) + 5 * (ip - 1) for it in 1:nT, ip in 1:nP]   # (T_idx, P)
    σf = Float32[(1 + 0.1f0 * kv) * (1 + 0.01f0 * it) * (1 + 0.001f0 * ip) *
                 Float32(exp(-((νf[i] - 4003.5) / 3)^2))
                 for i in 1:nν, kv in 1:nb, it in 1:nT, ip in 1:nP]               # (ν, b, T, P)

    fn = joinpath(mktempdir(), "co2_fixture.hdf")
    NCDataset(fn, "c") do ds
        defDim(ds, "nu", nν); defDim(ds, "b", nb); defDim(ds, "t", nT); defDim(ds, "p", nP)
        defDim(ds, "one", 1)
        defVar(ds, "Gas_Index", ["02"], ("one",))
        defVar(ds, "Wavenumber", νf, ("nu",))
        defVar(ds, "Pressure", pf .* 100, ("p",))                 # Pa
        defVar(ds, "Temperature", Tf, ("t", "p"))
        defVar(ds, "Broadener_01_VMR", vf, ("b",))
        defVar(ds, "Gas_02_Absorption", σf, ("nu", "b", "t", "p"))
    end

    lut = read_absco(fn)                                          # Float32, CPU
    @test lut isa AbscoLUT && lut.mol == 2
    @test lut.ν ≈ νf && lut.p ≈ pf && lut.T ≈ Tf && lut.vmr ≈ vf
    @test size(lut.σ) == (nν, nb, nT, nP)
    # Faithful: querying at every original node recovers the stored spectrum.
    for ip in 1:nP, it in 1:nT, kv in 1:nb
        @test compute_cross_section(lut, νf, pf[ip], Tf[it, ip]; vmr = vf[kv]) ≈ σf[:, kv, it, ip]
    end
end
