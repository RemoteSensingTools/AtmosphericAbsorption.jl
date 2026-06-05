using AtmosphericAbsorption
using DelimitedFiles: readdlm
using Test

# Validate HITRAN line-by-line cross-sections against HAPI golden spectra. Each golden
# ships the filtered .par subset (so we build the same line list) and HAPI's σ(ν) on a
# fixed grid; we compare in the line cores (the far wings underflow and are dominated by
# tiny absolute differences). Goldens are made by benchmark/hapi_reference/generate_golden.py
# and ship in the `reference_data` lazy artifact; GOLDEN (its golden dir) comes from runtests.jl.

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

# HAPI-style direct download from hitran.org. Opt-in (set AA_NETWORK_TESTS=1) so CI does
# not hit HITRAN's rate-limited API on every run; still skips gracefully if offline.
@testset "HITRAN direct download (network)" begin
    if !haskey(ENV, "AA_NETWORK_TESTS")
        @info "HITRAN download test skipped — set AA_NETWORK_TESTS=1 to enable"
    else
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

@testset "HITRAN extended params: line mixing + self-broadening (offline golden)" begin
    # Real CO2 4.3µm Q-branch slice with appended HT + Y + n_self/δ_self columns (carved
    # from a HITRANonline non-Voigt fetch). Verifies the parser fills Y_LM from
    # Y_HT_air_296 / Y_SDV_air_296 and n_self / δ_self from n_self / delta_self.
    parse_nv = AtmosphericAbsorption.LineLists._parse_nonvoigt_data
    path = joinpath(GOLDEN, "co2_linemix_2349_2351.data")
    db = parse_nv(path; numin = 2349.5, numax = 2351.0, min_strength = 0.0)
    @test count(!iszero, db.Y_LM) > 100          # most lines here carry Rosenkranz Y
    @test all(iszero, db.n_Y_LM)                 # HITRAN ships no Y temperature exponent
    @test count(db.n_self .!= db.n_air) > 100    # CO2 ships a distinct self T-exponent
    @test count(!iszero, db.δ_self) > 100        # …and a self pressure shift

    # Independent raw re-parse (Y_HT t8, Y_SDV t9, n_self t10, δ_self t11), ν-sorted.
    g(t, k) = (k ≤ length(t) && !(strip(t[k]) in ("", "#"))) ? parse(Float64, t[k]) : NaN
    raw = NamedTuple{(:ν, :Y, :ns, :ds),NTuple{4,Float64}}[]
    for ln in eachline(path)
        length(ln) < 160 && continue
        t = split(ln[161:end], ',')
        yht, ysdv = g(t, 8), g(t, 9)
        push!(raw, (ν = parse(Float64, strip(ln[4:15])),
                    Y = !isnan(yht) ? yht : (isnan(ysdv) ? 0.0 : ysdv),
                    ns = g(t, 10), ds = g(t, 11)))   # all present in this golden (no fallback)
    end
    sort!(raw, by = r -> r.ν)
    @test length(raw) == length(db)
    @test maximum(abs(raw[k].Y  - db.Y_LM[k])   for k in eachindex(db.Y_LM)) == 0.0
    @test maximum(abs(raw[k].ns - db.n_self[k]) for k in eachindex(db.n_self)) == 0.0
    @test maximum(abs(raw[k].ds - db.δ_self[k]) for k in eachindex(db.δ_self)) == 0.0

    # Line mixing must actually move the cross-section (asymmetric Q-branch redistribution).
    db0 = AtmosphericAbsorption.LineLists.LineDatabase(; mol = db.mol, iso = db.iso, ν0 = db.ν0,
        S = db.S, E_lower = db.E_lower, g_upper = db.g_upper, γ_air = db.γ_air, γ_self = db.γ_self,
        n_air = db.n_air, δ_air = db.δ_air, n_self = db.n_self, δ_self = db.δ_self,
        γ2_air = db.γ2_air, δ2_air = db.δ2_air, η = db.η,
        νVC = db.νVC, Y_LM = zero(db.Y_LM), molar_mass = db.molar_mass, meta = db.meta)
    grid = collect(2349.5:0.004:2351.0)
    mk(d) = LineByLineModel(d, TIPS2017PF(); profile = SpeedDependentVoigt(), wing_cutoff = 40.0)
    σmix = compute_cross_section(mk(db),  grid, 1013.25, 296.0)
    σno  = compute_cross_section(mk(db0), grid, 1013.25, 296.0)
    @test all(isfinite, σmix)
    @test maximum(abs.(σmix .- σno)) / maximum(σno) > 1e-3
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
