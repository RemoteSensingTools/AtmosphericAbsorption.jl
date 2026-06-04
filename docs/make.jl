using Documenter
using DocumenterVitepress
using AtmosphericAbsorption
using AtmosphericAbsorption.LineShapes: voigt, lorentz, doppler, pcqsdhc, HumlicekWeideman32

# ---------------------------------------------------------------------------
# Plots.jl with the plotly() backend → standalone interactive HTML per figure, embedded
# with an <iframe>. Real plots, computed from the package (data precomputed offline where a
# HITRAN key/network is needed, so the build stays CI-safe).
# ---------------------------------------------------------------------------
using Plots
plotly()
default(framestyle = :box, gridalpha = 0.25, linewidth = 2, size = (820, 460),
        fontfamily = "Inter, system-ui, sans-serif", legend = :topright)

const PLOTDIR = joinpath(@__DIR__, "src", "assets", "plots")
save_html(p, name) = (mkpath(PLOTDIR); savefig(p, joinpath(PLOTDIR, name * ".html")))
clip(y, fl = 1e-26) = max.(y, fl)   # keep log plots positive past first-order line-mixing wing negatives

# ---------------------------------------------------------------------------
# Generated figures
# ---------------------------------------------------------------------------
function plot_lineshape_families()
    cpf = HumlicekWeideman32()
    Δ = collect(-1.0:0.002:1.0)
    γd, γl = 0.08, 0.12
    y = sqrt(log(2.0)) * γl / γd
    p = plot(; title = "Area-normalized line-shape families (same γd, γl)",
             xlabel = "ν − ν₀  [cm⁻¹]", ylabel = "ϕ(ν)  [cm]")
    plot!(p, Δ, [doppler(d, γd) for d in Δ]; label = "Doppler")
    plot!(p, Δ, [lorentz(d, γl) for d in Δ]; label = "Lorentz")
    plot!(p, Δ, [voigt(cpf, d, γd, y) for d in Δ]; label = "Voigt")
    plot!(p, Δ, [real(pcqsdhc(cpf, 0.0, γd, γl, 0.03, 0.0, 0.0, 0.0, 0.0, d)) for d in Δ]; label = "Speed-dependent Voigt")
    plot!(p, Δ, [real(pcqsdhc(cpf, 0.0, γd, γl, 0.03, 0.0, 0.0, 0.05, 0.4, d)) for d in Δ]; label = "Hartmann-Tran")
    save_html(p, "lineshape_families")
end

# CO2 1.6 µm band at three temperatures — line strengths and widths shift with T.
function plot_temperature()
    par = joinpath(@__DIR__, "..", "test", "golden", "co2_6300_6400_p500_T250.par")
    db    = load_lines(HitranPort(par); mol = :CO2, iso = :main, min_strength = 1e-26)
    model = LineByLineModel(db; profile = Voigt(), wing_cutoff = 40.0)   # partition rides on db
    grid  = collect(6300.0:0.02:6400.0)
    p = plot(; title = "CO₂ 1.6 µm band vs temperature (p = 500 hPa)", yscale = :log10,
             xlabel = "wavenumber  [cm⁻¹]", ylabel = "σ  [cm²/molecule]", ylims = (1e-26, 1e-21))
    for T in (296.0, 250.0, 220.0)
        plot!(p, grid, clip(compute_cross_section(model, grid, 500.0, T)); label = "$(Int(T)) K")
    end
    save_html(p, "co2_temperature")
end

# Read a committed cross-section overlay (cols: ν, σ_a, σ_b). Generated offline by the
# benchmark/gen_*.jl scripts so the doc build needs no network/HITRAN key (CI-safe).
function read_xsec(name)
    nu, a, b = Float64[], Float64[], Float64[]
    for ln in eachline(joinpath(@__DIR__, "src", "assets", name))
        startswith(strip(ln), "#") && continue
        c = split(ln)
        isempty(c) && continue
        push!(nu, parse(Float64, c[1])); push!(a, parse(Float64, c[2])); push!(b, parse(Float64, c[3]))
    end
    nu, a, b
end

# ExoMol vs HITRAN CO cross-section — ExoMol derives line strengths from Einstein-A
# coefficients + its own partition function, yet lands on top of HITRAN.
function plot_exomol_co()
    ν, σe, σh = read_xsec("exomol_co_xsec.txt")
    p = plot(; title = "CO cross-section — ExoMol vs HITRAN (p = 1013 hPa, T = 296 K)", yscale = :log10,
             xlabel = "wavenumber  [cm⁻¹]", ylabel = "σ  [cm²/molecule]")
    plot!(p, ν, clip(σe); label = "ExoMol Li2015 (S from Einstein-A)")
    plot!(p, ν, clip(σh); label = "HITRAN", linestyle = :dot)
    save_html(p, "exomol_co_xsec")
end

# CO2 4.3 µm band: first-order line mixing redistributes intensity between overlapping
# lines and suppresses the far wings (sub-Lorentzian) vs a plain Voigt sum.
function plot_co2_linemix()
    ν, σno, σlm = read_xsec("co2_linemix_band.txt")
    p = plot(; title = "CO₂ 4.3 µm band — effect of line mixing (p = 1013 hPa, T = 296 K)", yscale = :log10,
             xlabel = "wavenumber  [cm⁻¹]", ylabel = "σ  [cm²/molecule]")
    plot!(p, ν, clip(σno); label = "Voigt (no line mixing)")
    plot!(p, ν, clip(σlm); label = "with line mixing (Hartmann-Tran)")
    save_html(p, "co2_linemix")
end

# H2O: the speed dependence of collisions (Hartmann-Tran γ₂/ν_VC) narrows and reshapes the
# line cores relative to a plain Voigt.
function plot_h2o_ht()
    ν, σv, σht = read_xsec("h2o_voigt_vs_ht.txt")
    p = plot(; title = "H₂O speed dependence — Voigt vs Hartmann-Tran (p = 1013 hPa, T = 250 K)", yscale = :log10,
             xlabel = "wavenumber  [cm⁻¹]", ylabel = "σ  [cm²/molecule]")
    plot!(p, ν, clip(σv);  label = "Voigt")
    plot!(p, ν, clip(σht); label = "Hartmann-Tran (γ₂, ν_VC)", linestyle = :dash)
    save_html(p, "h2o_voigt_vs_ht")
end

function plot_benchmark()
    cut = ["2.5", "5", "10", "25", "full"]
    x = 1:5
    p = plot(; title = "Time vs wing cutoff — 4000 lines × 400 cm⁻¹ band (A100)", yscale = :log10,
             xlabel = "wing cutoff  [cm⁻¹]", ylabel = "time  [ms]", xticks = (x, cut), marker = :circle)
    plot!(p, x, [114.0, 173.0, 286.0, 623.0, 4691.0]; label = "hapi2 (numba)")
    plot!(p, x, [177.0, 235.0, 349.0, 683.0, 4122.0]; label = "AtmosphericAbsorption CPU")
    plot!(p, x, [2.1, 2.2, 2.4, 3.0, 10.1]; label = "AtmosphericAbsorption GPU")
    save_html(p, "benchmark")
end

# Generate the species/isotopologue reference table from the bundled data + registry, so
# the docs page stays in sync with the code that resolves :CO2 / :main / formulas.
function write_isotopologue_table()
    csv = joinpath(@__DIR__, "..", "data", "hitran_isotopologues.csv")
    rows = sort!([(mol = parse(Int, a[1]), iso = parse(Int, a[2]), formula = a[3])
                  for ln in eachline(csv) if !startswith(ln, "mol")
                  for a in (split(strip(ln), ','),) if length(a) ≥ 3], by = r -> (r.mol, r.iso))
    fmt(f) = try string(round(f(); sigdigits = 5)) catch; "—" end
    io = IOBuffer()
    println(io, "# Species & isotopologues\n")
    println(io, "Every molecule below can be named with the [species notation](data_sources.md#species-notation) — ",
                "`mol = :CO2` / `\"CO2\"` / `2`, and isotopologues with `iso = :main` / an id / `:ALL`. ",
                "The mapping is overridable with `register_molecule!` / `register_isotopologue!`. ",
                "Molar mass / abundance are blank for the heavy molecules that ship only as ",
                "[tabulated cross-sections](data_sources.md).\n")
    println(io, "| id | molecule | iso | isotopologue | molar mass [g/mol] | abundance |")
    println(io, "|---:|:---------|----:|:-------------|-------------------:|----------:|")
    for r in rows
        sym = string(AtmosphericAbsorption.molecule_symbol(r.mol))
        mm  = fmt(() -> AtmosphericAbsorption.LineLists.molar_mass(r.mol, r.iso))
        ab  = fmt(() -> AtmosphericAbsorption.LineLists.abundance(r.mol, r.iso))
        println(io, "| $(r.mol) | `:$sym` | $(r.iso) | `$(r.formula)` | $mm | $ab |")
    end
    write(joinpath(@__DIR__, "src", "isotopologues.md"), String(take!(io)))
end

@info "Generating documentation figures…"
plot_lineshape_families()
plot_temperature()
plot_co2_linemix()
plot_h2o_ht()
plot_exomol_co()
plot_benchmark()
write_isotopologue_table()

const REPO = "github.com/RemoteSensingTools/AtmosphericAbsorption.jl"
const IN_CI = get(ENV, "CI", "false") == "true"

makedocs(;
    sitename = "AtmosphericAbsorption.jl",
    modules = [AtmosphericAbsorption],
    clean = true,
    checkdocs = :none,
    format = MarkdownVitepress(;
        repo = REPO,
        devbranch = "main",
        devurl = "dev",
        # Deploy only in CI; local builds just render.
        deploy_decision = IN_CI ? nothing : Documenter.DeployDecision(all_ok = false),
        description = "GPU-accelerated molecular absorption cross-sections for atmospheric radiative transfer.",
        assets = ["assets/logo.png"],
        sidebar_drawer = true,
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Line shapes" => "line_shapes.md",
        "Data sources" => "data_sources.md",
        "Species & isotopologues" => "isotopologues.md",
        "Continuum" => "continuum.md",
        "GPU & precision" => "gpu.md",
        "Benchmarks" => "benchmarks.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
)

if IN_CI
    DocumenterVitepress.deploydocs(;
        root = @__DIR__,
        repo = "$(REPO).git",
        target = "build",
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true,
    )
end
