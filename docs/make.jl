using Documenter
using DocumenterVitepress
using AtmosphericAbsorption
using AtmosphericAbsorption.LineShapes: voigt, lorentz, doppler, pcqsdhc, HumlicekWeideman32

# ---------------------------------------------------------------------------
# Minimal JSON serialization + Plotly asset writer (adapted from vSmartMOM's docs).
# Each plot is a standalone HTML file under assets/plots/, embedded with an <iframe>.
# ---------------------------------------------------------------------------
const PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"

_jesc(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
_json(x::Nothing) = "null"
_json(x::Bool) = x ? "true" : "false"
_json(x::Integer) = string(x)
_json(x::AbstractFloat) = isfinite(x) ? string(x) : "null"
_json(x::AbstractString) = "\"" * _jesc(x) * "\""
_json(x::Symbol) = _json(String(x))
_json(x::AbstractVector) = "[" * join(_json.(x), ",") * "]"
_json(x::Tuple) = "[" * join(_json.(collect(x)), ",") * "]"
_json(x::NamedTuple) = "{" * join((_json(String(k)) * ":" * _json(v) for (k, v) in pairs(x)), ",") * "}"
_json(x::AbstractDict) = "{" * join((_json(String(k)) * ":" * _json(v) for (k, v) in x), ",") * "}"

function write_plot(name, data, layout)
    dir = joinpath(@__DIR__, "src", "assets", "plots")
    mkpath(dir)
    config = (; responsive = true, displaylogo = false)
    html = """
    <!doctype html><html><head><meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <script src="$PLOTLY_CDN"></script>
    <style>html,body,#p{width:100%;height:100%;margin:0;background:#fff}</style></head>
    <body><div id="p"></div><script>
    Plotly.newPlot("p", $(_json(data)), $(_json(layout)), $(_json(config)));
    </script></body></html>
    """
    write(joinpath(dir, name * ".html"), html)
end

trace(x, y, nm) = (; x = collect(Float64, x), y = collect(Float64, y), mode = "lines",
                   name = nm, type = "scatter")
axis(t) = (; title = t, showgrid = true, gridcolor = "#eee", zeroline = false)
lay(title, xt, yt; ytype = "linear") =
    (; title = (; text = title), xaxis = axis(xt), yaxis = merge(axis(yt), (; type = ytype)),
     legend = (; orientation = "h", y = -0.2), margin = (; l = 70, r = 20, t = 50, b = 60),
     font = (; family = "Inter, system-ui, sans-serif"), paper_bgcolor = "#fff", plot_bgcolor = "#fff")

# ---------------------------------------------------------------------------
# Generated figures
# ---------------------------------------------------------------------------
function plot_lineshape_families()
    cpf = HumlicekWeideman32()
    Δ = collect(-1.0:0.002:1.0)
    γd, γl = 0.08, 0.12
    y = sqrt(log(2.0)) * γl / γd
    vg = [voigt(cpf, d, γd, y) for d in Δ]
    lo = [lorentz(d, γl) for d in Δ]
    dp = [doppler(d, γd) for d in Δ]
    # speed-dependent Voigt (Γ2≠0) and Hartmann-Tran (νVC, η ≠ 0) via pcqsdhc
    sdv = [real(pcqsdhc(cpf, 0.0, γd, γl, 0.03, 0.0, 0.0, 0.0, 0.0, d)) for d in Δ]
    ht  = [real(pcqsdhc(cpf, 0.0, γd, γl, 0.03, 0.0, 0.0, 0.05, 0.4, d)) for d in Δ]
    write_plot("lineshape_families",
        [trace(Δ, dp, "Doppler"), trace(Δ, lo, "Lorentz"), trace(Δ, vg, "Voigt"),
         trace(Δ, sdv, "Speed-dependent Voigt"), trace(Δ, ht, "Hartmann-Tran")],
        lay("Area-normalized line-shape families (same γd, γl)", "ν − ν₀  [cm⁻¹]", "ϕ(ν)  [cm]"))
end

function plot_co2_crosssection()
    par = joinpath(@__DIR__, "..", "test", "golden", "co2_6300_6400_p500_T250.par")
    db = load_lines(HitranPort(par); mol = :CO2, iso = :main, min_strength = 1e-26)
    model = LineByLineModel(db; profile = Voigt(), wing_cutoff = 40.0)   # partition rides on db
    grid = collect(6300.0:0.01:6400.0)
    σ250 = compute_cross_section(model, grid, 500.0, 250.0)
    σ296 = compute_cross_section(model, grid, 500.0, 296.0)
    write_plot("co2_crosssection",
        [trace(grid, σ296, "296 K"), trace(grid, σ250, "250 K")],
        lay("CO₂ absorption cross-section, 1.6 µm band (p = 500 hPa)",
            "wavenumber  [cm⁻¹]", "σ  [cm²/molecule]"; ytype = "log"))
end

# Read a committed cross-section overlay (cols: ν, σ_a, σ_b). Generated offline by
# benchmark/gen_exomol_fig.jl so the doc build needs no network (CI-safe).
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
    nu, σe, σh = read_xsec("exomol_co_xsec.txt")
    te = (; x = nu, y = σe, mode = "lines", type = "scatter",
          name = "ExoMol Li2015 (S from Einstein-A)", line = (; color = "#1f77b4", width = 2))
    th = (; x = nu, y = σh, mode = "lines", type = "scatter",
          name = "HITRAN", line = (; color = "#d62728", width = 2, dash = "dot"))
    write_plot("exomol_co_xsec", [te, th],
        lay("CO cross-section — ExoMol vs HITRAN (p = 1013 hPa, T = 296 K)",
            "wavenumber  [cm⁻¹]", "σ  [cm²/molecule]"; ytype = "log"))
end

function plot_benchmark()
    cut = ["2.5", "5", "10", "25", "full"]
    hapi2 = [114.0, 173.0, 286.0, 623.0, 4691.0]
    cpu   = [177.0, 235.0, 349.0, 683.0, 4122.0]
    gpu   = [2.1, 2.2, 2.4, 3.0, 10.1]
    bar(y, nm) = (; x = cut, y = y, type = "bar", name = nm)
    write_plot("benchmark",
        [bar(hapi2, "hapi2 (numba)"), bar(cpu, "AtmosphericAbsorption CPU"),
         bar(gpu, "AtmosphericAbsorption GPU")],
        (; title = (; text = "Time vs wing cutoff — 4000 lines × 400 cm⁻¹ band (A100)"),
         barmode = "group", xaxis = axis("wing cutoff  [cm⁻¹]"),
         yaxis = merge(axis("time  [ms]"), (; type = "log")),
         legend = (; orientation = "h", y = -0.2), margin = (; l = 70, r = 20, t = 50, b = 60),
         font = (; family = "Inter, system-ui, sans-serif"), paper_bgcolor = "#fff", plot_bgcolor = "#fff"))
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
plot_co2_crosssection()
write_isotopologue_table()
plot_exomol_co()
plot_benchmark()

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
