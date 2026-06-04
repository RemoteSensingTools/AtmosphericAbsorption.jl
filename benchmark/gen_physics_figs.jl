# Precompute the advanced-physics figure data for the docs (committed so the build stays
# offline): CO2 4.3µm band with/without line mixing, and H2O Voigt vs Hartmann-Tran.
# Run: source ~/.bashrc; julia --project=. benchmark/gen_physics_figs.jl   (needs HITRAN key)
using AtmosphericAbsorption
const LL = AtmosphericAbsorption.LineLists
activate_hitran!(get(ENV, "HITRAN_API_KEY", ""))
const ASSET = joinpath(@__DIR__, "..", "docs", "src", "assets")

writecols(name, hdr, ν, cols...) = open(joinpath(ASSET, name), "w") do io
    println(io, "# ", hdr)
    for i in eachindex(ν)
        print(io, ν[i]); for c in cols; print(io, "  ", c[i]); end; println(io)
    end
end

# --- CO2 4.3 µm band: line mixing on vs off (HITRAN gives SDV line-mixing Y; no γ2) ---
@info "CO2 4.3µm…"
co2 = load_hitran_nonvoigt(:CO2; numin = 2280, numax = 2400, min_strength = 1e-25)
@info "CO2 lines" n = length(co2) nY = count(!iszero, co2.Y_LM) nγ2 = count(!iszero, co2.γ2_air)
co2_nomix = LL.LineDatabase(; mol = co2.mol, iso = co2.iso, ν0 = co2.ν0, S = co2.S,
    E_lower = co2.E_lower, g_upper = co2.g_upper, γ_air = co2.γ_air, γ_self = co2.γ_self,
    n_air = co2.n_air, δ_air = co2.δ_air, n_self = co2.n_self, δ_self = co2.δ_self,
    γ2_air = co2.γ2_air, δ2_air = co2.δ2_air, νVC = co2.νVC, η = co2.η, n_γ2 = co2.n_γ2,
    Y_LM = zero(co2.Y_LM), molar_mass = co2.molar_mass, meta = co2.meta, partition = co2.partition)
g = collect(2280.0:0.02:2400.0)
σmix = compute_cross_section(LineByLineModel(co2;       profile = HartmannTran(), wing_cutoff = 60.0), g, 1013.25, 296.0)
σno  = compute_cross_section(LineByLineModel(co2_nomix; profile = Voigt(),        wing_cutoff = 60.0), g, 1013.25, 296.0)
writecols("co2_linemix_band.txt", "CO2 4.3µm, p=1013hPa T=296K. nu  sigma_Voigt_noLM  sigma_lineMixing", g, σno, σmix)

# --- H2O speed dependence: Voigt vs Hartmann-Tran (same lines, HT uses γ2/νVC) ---
@info "H2O HT…"
h2o = load_hitran_nonvoigt(:H2O; numin = 3700, numax = 3850, min_strength = 1e-25)
@info "H2O lines" n = length(h2o) nγ2 = count(!iszero, h2o.γ2_air)
gh = collect(3800.0:0.004:3825.0)
σv  = compute_cross_section(LineByLineModel(h2o; profile = Voigt(),        wing_cutoff = 40.0), gh, 1013.25, 250.0)
σht = compute_cross_section(LineByLineModel(h2o; profile = HartmannTran(), wing_cutoff = 40.0), gh, 1013.25, 250.0)
writecols("h2o_voigt_vs_ht.txt", "H2O 3800-3825, p=1013hPa T=250K. nu  sigma_Voigt  sigma_HartmannTran", gh, σv, σht)
println("DONE  CO2 npts=", length(g), "  H2O npts=", length(gh))
