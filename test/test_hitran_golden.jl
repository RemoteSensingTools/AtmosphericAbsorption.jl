using AtmosphericAbsorption
using DelimitedFiles: readdlm
using Test

# Validate HITRAN line-by-line cross-sections against HAPI golden spectra. Each golden
# ships the filtered .par subset (so we build the same line list) and HAPI's σ(ν) on a
# fixed grid; we compare in the line cores (the far wings underflow and are dominated by
# tiny absolute differences). Goldens are made by benchmark/hapi_reference/generate_golden.py.
const GOLDEN = joinpath(@__DIR__, "golden")

# name, mol, iso, p[hPa], T[K], wing[cm⁻¹]
const GOLDEN_CASES = [
    ("co2_6300_6400_p500_T250", 2, 1, 500.0, 250.0, 40.0),
    ("co2_6300_6400_p1013_T296", 2, 1, 1013.25, 296.0, 40.0),
    ("h2o_7000_7100_p800_T280", 1, 1, 800.0, 280.0, 40.0),
]

@testset "HITRAN local_iso_id decoding" begin
    f = AtmosphericAbsorption.LineLists._isoid
    @test (f('1'), f('9'), f('0'), f('A'), f('B')) == (1, 9, 10, 11, 12)
end

@testset "HITRAN API key management" begin
    activate_hitran!("test-dummy-key")          # never a real key in tests
    @test AtmosphericAbsorption.LineLists.has_hitran_api_key()
    @test AtmosphericAbsorption.LineLists.hitran_api_key() == "test-dummy-key"
    AtmosphericAbsorption.LineLists._HITRAN_API_KEY[] = nothing   # reset session key
end

# HAPI-style direct download from hitran.org — network-gated (skipped, not failed, offline).
@testset "HITRAN direct download (network)" begin
    db = try
        load_lines(HitranPort(; molecule = "CO", numin = 2100.0, numax = 2200.0);
                   mol = 5, iso = 1, min_strength = 1e-25)
    catch e
        @info "HITRAN download failed (offline?) — skipping" exception = (e, catch_backtrace())
        nothing
    end
    if db !== nothing
        @test length(db) > 40
        @test all(2100 .≤ db.ν0 .≤ 2200)
        @test issorted(db.ν0)
    end
end

# Authenticated non-Voigt (HT/SDV) fetch — gated on BOTH an API key and network.
@testset "HITRAN non-Voigt (HT/SDV) fetch" begin
    if !AtmosphericAbsorption.LineLists.has_hitran_api_key()
        @info "No HITRAN API key (activate_hitran!/HITRAN_API_KEY) — skipping non-Voigt test"
    else
        db = try
            load_hitran_nonvoigt("H2O"; numin = 3700.0, numax = 3850.0, min_strength = 1e-25)
        catch e
            @info "HITRAN non-Voigt fetch failed — skipping" exception = (e, catch_backtrace())
            nothing
        end
        if db !== nothing
            @test count(!iszero, db.γ2_air) > 100        # ~336 H2O HT lines in this band
            grid = collect(3700.0:0.05:3850.0)
            σsdv = compute_cross_section(LineByLineModel(db, TIPS2017PF(); profile = SpeedDependentVoigt(), wing_cutoff = 40.0), grid, 500.0, 250.0)
            σvgt = compute_cross_section(LineByLineModel(db, TIPS2017PF(); profile = Voigt(), wing_cutoff = 40.0), grid, 500.0, 250.0)
            @test all(isfinite, σsdv) && maximum(σsdv) > 0
            @test maximum(abs.(σsdv .- σvgt)) / maximum(σvgt) > 1e-3   # HT params change the shape
        end
    end
end

@testset "HITRAN cross-section vs HAPI golden" begin
    @testset "$name" for (name, mol, iso, p, T, wing) in GOLDEN_CASES
        d = readdlm(joinpath(GOLDEN, name * ".txt"), comments = true, comment_char = '#')
        ν, σ_ref = d[:, 1], d[:, 2]
        db = load_lines(HitranPort(joinpath(GOLDEN, name * ".par")); mol = mol, iso = iso)
        model = LineByLineModel(db, TIPS2017PF(); profile = Voigt(), wing_cutoff = wing)
        σ = compute_cross_section(model, collect(ν), p, T)
        mask = σ_ref .> 1e-3 * maximum(σ_ref)          # compare in the line cores
        core_relerr = maximum(abs.(σ[mask] .- σ_ref[mask]) ./ σ_ref[mask])
        @test core_relerr < 5e-3
    end
end
