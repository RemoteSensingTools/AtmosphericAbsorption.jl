using AtmosphericAbsorption
using AtmosphericAbsorption.LineLists: parse_def, read_states, read_pf, read_broad,
                                       _broad_key, exomol_intensity, _stream_trans!
using Test

const FIX = joinpath(@__DIR__, "fixtures", "exomol")
const COPAR = joinpath(@__DIR__, "golden", "co_2100_2200.par")

@testset "ExoMol readers (synthetic fixtures)" begin
    def = parse_def(joinpath(FIX, "toy.def"))
    @test def.mass == 28.0
    @test def.nstates == 4
    @test (def.default_γ, def.default_n) == (0.07, 0.5)

    st = read_states(joinpath(FIX, "toy.states"), 4)
    @test st.E[4] == 2010.0 && st.g[4] == 3 && st.J[2] == 1

    pf = read_pf(joinpath(FIX, "toy.pf"))
    @test pf.T == [100.0, 296.0, 400.0] && pf.Q == [5.0, 10.0, 12.0]

    diet, tbl = read_broad(joinpath(FIX, "toy__air.broad"))
    @test diet == "m0"
    @test tbl[1] == (0.05, 0.7) && tbl[-1] == (0.06, 0.6)
    @test read_broad(joinpath(FIX, "does_not_exist.broad")) == ("", Dict{Int,Tuple{Float64,Float64}}())
end

@testset "ExoMol branch index (diet codes)" begin
    @test _broad_key("m0", 0, 1) == 1     # R branch (ΔJ=+1): m = J_lo+1
    @test _broad_key("m0", 5, 4) == -5    # P branch (ΔJ=-1): m = -J_lo
    @test _broad_key("m0", 3, 3) == 3     # Q branch (ΔJ=0):  m = J_lo
    @test _broad_key("a0", 7, 8) == 7     # a0 keys on lower J directly
end

@testset "ExoMol .trans streaming (offline fixture)" begin
    def = parse_def(joinpath(FIX, "toy.def"))
    st  = read_states(joinpath(FIX, "toy.states"), def.nstates)
    diet, broad = read_broad(joinpath(FIX, "toy__air.broad"))
    abund, Q296 = 0.5, 10.0     # arbitrary; toy.pf has Q(296)=10
    cols = (ν0 = Float64[], S = Float64[], E = Float64[], g = Float64[], γ = Float64[], n = Float64[])
    open(joinpath(FIX, "toy.trans")) do io
        _stream_trans!(cols, io, st, diet, broad, def, abund, Q296, 1000.0, 3000.0, 0.0)
    end
    @test cols.ν0 == [2010.0, 1990.0]                         # ν=10 line excluded by window
    @test cols.S[1] ≈ exomol_intensity(3, 5.0, 2010.0, 0.0, Q296, abund)
    @test cols.γ[1] ≈ 0.05 * 1.01325 && cols.n[1] == 0.7      # R branch (m=+1) bar→atm
    @test cols.n[2] == 0.6                                    # P branch (m=-1)
end

@testset "ExoMol intensity formula vs HITRAN CO golden" begin
    # CO 1-0 fundamental lines (g_up, A, ν, E_lo, Q296=107.4198, abund=0.986544) → HITRAN S.
    @test exomol_intensity(3, 11.95, 2147.0811, 0.0, 107.4198, 0.986544) ≈ 9.480e-20 rtol = 2e-3
    @test exomol_intensity(17, 17.52, 2172.7588, 107.642, 107.4198, 0.986544) ≈ 4.556e-19 rtol = 2e-3
end

# Full integration: fetch real ExoMol CO and validate against HITRAN. Network-gated —
# skipped (not failed) when exomol.com is unreachable.
@testset "ExoMol CO vs HITRAN (network)" begin
    port = ExoMolPort("CO", "12C-16O", "Li2015"; hitran_mol = 5, hitran_iso = 1)
    edb = try
        load_lines(port; ν_min = 2100.0, ν_max = 2200.0, min_strength = 1e-25)
    catch e
        @info "ExoMol fetch failed (offline?) — skipping network test" exception = (e, catch_backtrace())
        nothing
    end
    if edb !== nothing
        hdb = load_lines(HitranPort(COPAR); mol = 5, iso = 1, ν_min = 2100.0, ν_max = 2200.0, min_strength = 1e-25)
        @test length(edb) == length(hdb)
        maxdν, maxrelS = 0.0, 0.0
        for k in eachindex(edb.ν0)
            j = argmin(abs.(hdb.ν0 .- edb.ν0[k]))
            abs(hdb.ν0[j] - edb.ν0[k]) > 0.02 && continue
            maxdν = max(maxdν, abs(hdb.ν0[j] - edb.ν0[k]))
            maxrelS = max(maxrelS, abs(edb.S[k] - hdb.S[j]) / hdb.S[j])
        end
        @test maxdν < 0.01            # line positions
        @test maxrelS < 0.01          # intensities (the verified-formula deliverable)
    end
end
