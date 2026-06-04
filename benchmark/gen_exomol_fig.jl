# Regenerate the ExoMol-vs-HITRAN CO cross-section overlay used in the docs
# (docs/src/assets/exomol_co_xsec.txt). Fetches ExoMol CO once (network), loads the
# bundled HITRAN CO golden, and writes both cross-sections on a shared grid so the
# documentation build stays offline. Run: julia --project=. benchmark/gen_exomol_fig.jl
using AtmosphericAbsorption

const PKG = dirname(@__DIR__)

@info "Fetching ExoMol CO Li2015 (2100–2200 cm⁻¹)…"
eport = ExoMolPort("CO", "12C-16O", "Li2015"; hitran_mol = 5, hitran_iso = 1)
edb   = load_lines(eport; ν_min = 2100.0, ν_max = 2200.0, min_strength = 1e-25)
epf   = partition_function(eport, 5, 1)

@info "Loading HITRAN CO golden…"
hport = HitranPort(joinpath(PKG, "test", "golden", "co_2100_2200.par"))
hdb   = load_lines(hport; mol = 5, iso = 1, ν_min = 2100.0, ν_max = 2200.0, min_strength = 1e-25)
hpf   = partition_function(hport, 5, 1)

# matched line-intensity comparison (the ExoMol Einstein-A → S validation)
maxrel, nmatch = 0.0, 0
for k in eachindex(edb.ν0)
    j = argmin(abs.(hdb.ν0 .- edb.ν0[k]))
    abs(hdb.ν0[j] - edb.ν0[k]) > 0.02 && continue
    global nmatch += 1
    global maxrel = max(maxrel, abs(edb.S[k] - hdb.S[j]) / hdb.S[j])
end
@info "Matched lines" n_exomol=length(edb) n_hitran=length(hdb) nmatch maxrel_pct=100maxrel

# cross-section overlay on a shared sub-band, p=1013.25 hPa, T=296 K
emodel = LineByLineModel(edb, epf; profile = Voigt(), wing_cutoff = 40.0)
hmodel = LineByLineModel(hdb, hpf; profile = Voigt(), wing_cutoff = 40.0)
grid = collect(2140.0:0.01:2160.0)
σe = compute_cross_section(emodel, grid, 1013.25, 296.0)
σh = compute_cross_section(hmodel, grid, 1013.25, 296.0)

out = joinpath(PKG, "docs", "src", "assets", "exomol_co_xsec.txt")
open(out, "w") do io
    println(io, "# CO cross-section, p=1013.25 hPa, T=296 K. cols: nu[cm-1]  sigma_ExoMol  sigma_HITRAN  [cm^2/molecule]")
    println(io, "# ExoMol Li2015 (S from Einstein-A) vs HITRAN; matched-line max rel intensity error = ",
            round(100maxrel, sigdigits = 3), "% over ", nmatch, " lines")
    for i in eachindex(grid)
        println(io, grid[i], "  ", σe[i], "  ", σh[i])
    end
end
@info "Wrote" out npts=length(grid)
