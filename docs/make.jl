using Documenter
using DocumenterVitepress
using AtmosphericAbsorption
using AtmosphericAbsorption.LineShapes: voigt, lorentz, doppler, pcqsdhc, HumlicekWeideman32
using LazyArtifacts
import Pkg.Artifacts: ensure_artifact_installed

# Precomputed figure data + the golden CO2 .par the temperature plot uses live in the
# `reference_data` lazy artifact (out of the package tree). Resolve via the repo-root
# Artifacts.toml explicitly — `@artifact_str` would walk up from docs/make.jl, hit
# docs/Project.toml, and stop before reaching it. Downloads on demand; CI has network.
const REFDATA = try
    ensure_artifact_installed("reference_data", joinpath(@__DIR__, "..", "Artifacts.toml"))
catch err
    @warn "reference_data artifact unavailable — data-backed figures skipped" exception = err
    nothing
end

# ---------------------------------------------------------------------------
# Interactive figures: one self-contained Plotly HTML per figure, embedded in the pages with
# an <iframe>. We hand-write the HTML (CDN Plotly + a responsive, container-filling div) and
# drop it straight into the BUILT site (see `_asset_destinations`), because DocumenterVitepress
# only promotes logo/favicon from `assets/` into Vitepress's `public/` — any other static file
# is dropped. This mirrors vSmartMOM's docs. Figures are written AFTER `makedocs()`.
# ---------------------------------------------------------------------------
const PLOTLY_CDN_URL = "https://cdn.plot.ly/plotly-2.35.2.min.js"
const PAL = ("#2563eb", "#dc2626", "#16a34a", "#7c3aed", "#ea580c", "#0891b2")
const SIG = "σ  [cm²/molecule]"
const WAVENO = "wavenumber  [cm⁻¹]"
clip(y, fl = 1e-26) = max.(y, fl)   # keep log plots positive past first-order line-mixing wing negatives

_json_escape(s::AbstractString) =
    replace(s, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")

_json_value(::Nothing) = "null"
_json_value(x::Bool) = x ? "true" : "false"
_json_value(x::Integer) = string(x)
_json_value(x::AbstractFloat) = isfinite(x) ? string(x) : "null"
_json_value(x::AbstractString) = "\"" * _json_escape(x) * "\""
_json_value(x::Symbol) = _json_value(String(x))
_json_value(x::Union{Tuple,AbstractVector}) = "[" * join((_json_value(v) for v in x), ",") * "]"
_json_value(x::NamedTuple) = "{" * join((_json_value(String(k)) * ":" * _json_value(v) for (k, v) in pairs(x)), ",") * "}"

# Every place a built copy of `assets/<subdir>/` must land for the docs to serve it: the
# Vitepress source tree, its `public/`, and each built base (read from `bases.txt`). Valid only
# after `makedocs()` has produced `build/`.
function _asset_destinations(subdir)
    dests = [joinpath(@__DIR__, "build", ".documenter", "assets", subdir),
             joinpath(@__DIR__, "build", ".documenter", "public", "assets", subdir)]
    bases = joinpath(@__DIR__, "build", "bases.txt")
    isfile(bases) && append!(dests,
        [joinpath(@__DIR__, "build", string(i), "assets", subdir) for i in eachindex(readlines(bases))])
    return unique(dests)
end

function _plotly_doc(data, layout; title)
    config = (; responsive = true, displaylogo = false, modeBarButtonsToRemove = ["lasso2d", "select2d"])
    return """
    <!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <title>$title</title>
      <script src="$PLOTLY_CDN_URL"></script>
      <style>html, body, #plot { width: 100%; height: 100%; margin: 0; background: transparent; }</style>
    </head>
    <body>
      <div id="plot"></div>
      <script>
        // Transparent background lets the embedding docs page show through, so the figure blends
        // into both light and dark mode; we adapt text / gridlines / axes for readability. The
        // figure is same-origin with the page, so it reads VitePress's `html.dark` and follows the
        // day/night toggle live (falling back to prefers-color-scheme if ever cross-origin).
        const data = $(_json_value(data)), config = $(_json_value(config)), base = $(_json_value(layout));
        const THEME = {
          light: { font: "#3c3c43", grid: "rgba(60,60,67,0.14)", zero: "rgba(60,60,67,0.30)", axis: "rgba(60,60,67,0.30)" },
          dark:  { font: "#dfdfd6", grid: "rgba(235,235,245,0.12)", zero: "rgba(235,235,245,0.26)", axis: "rgba(235,235,245,0.22)" }
        };
        function isDark() {
          try { return !!window.parent.document.documentElement.classList.contains("dark"); }
          catch (e) { return !!(window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches); }
        }
        function attrs(t) {
          return {
            paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
            "font.color": t.font, "legend.font.color": t.font,
            "xaxis.gridcolor": t.grid, "yaxis.gridcolor": t.grid,
            "xaxis.zerolinecolor": t.zero, "yaxis.zerolinecolor": t.zero,
            "xaxis.linecolor": t.axis, "yaxis.linecolor": t.axis,
            "xaxis.tickcolor": t.axis, "yaxis.tickcolor": t.axis
          };
        }
        // Bake the initial theme into the layout so there is no first-paint flash, then relayout on toggle.
        const t0 = isDark() ? THEME.dark : THEME.light;
        const layout = Object.assign({}, base, {
          paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
          font: Object.assign({}, base.font || {}, { color: t0.font }),
          xaxis: Object.assign({}, base.xaxis || {}, { gridcolor: t0.grid, zerolinecolor: t0.zero, linecolor: t0.axis, tickcolor: t0.axis }),
          yaxis: Object.assign({}, base.yaxis || {}, { gridcolor: t0.grid, zerolinecolor: t0.zero, linecolor: t0.axis, tickcolor: t0.axis }),
          legend: Object.assign({}, base.legend || {}, { font: { color: t0.font } })
        });
        Plotly.newPlot("plot", data, layout, config);
        function sync() { Plotly.relayout("plot", attrs(isDark() ? THEME.dark : THEME.light)); }
        try { new MutationObserver(sync).observe(window.parent.document.documentElement, { attributes: true, attributeFilter: ["class"] }); }
        catch (e) { if (window.matchMedia) matchMedia("(prefers-color-scheme: dark)").addEventListener("change", sync); }
      </script>
    </body>
    </html>
    """
end

function _write_plotly_asset(filename, data, layout; title)
    html = _plotly_doc(data, layout; title)
    for dest in _asset_destinations("plots")
        mkpath(dest)
        write(joinpath(dest, filename), html)
    end
    return nothing
end

_standard_font() = (; family = "Inter, system-ui, sans-serif")
_bottom_legend() = (; orientation = "h", x = 0.0, y = -0.22, xanchor = "left", yanchor = "top")

function _line(name, x, y; color = nothing, dash = nothing)
    line = (; width = 2.5)
    isnothing(color) || (line = merge(line, (; color)))
    isnothing(dash) || (line = merge(line, (; dash)))
    return (; type = "scatter", mode = "lines", name, x = collect(Float64, x), y = collect(Float64, y), line)
end

function _layout(title, xtitle, ytitle; ylog = false, xcategory = false)
    return (; title = (; text = title, x = 0.02),
            margin = (; l = 74, r = 20, t = 70, b = 86),
            xaxis = (; title = xtitle, type = xcategory ? "category" : "-"),
            yaxis = (; title = ytitle, type = ylog ? "log" : "linear"),
            legend = _bottom_legend(),
            paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", font = _standard_font())
end

# ---------------------------------------------------------------------------
# Generated figures
# ---------------------------------------------------------------------------
function plot_lineshape_families()
    cpf = HumlicekWeideman32()
    Δ = collect(-1.0:0.002:1.0)
    γd, γl = 0.08, 0.12
    y = sqrt(log(2.0)) * γl / γd
    data = [
        _line("Doppler", Δ, [doppler(d, γd) for d in Δ]; color = PAL[1]),
        _line("Lorentz", Δ, [lorentz(d, γl) for d in Δ]; color = PAL[2]),
        _line("Voigt", Δ, [voigt(cpf, d, γd, y) for d in Δ]; color = PAL[3]),
        _line("Speed-dependent Voigt", Δ,
              [real(pcqsdhc(cpf, 0.0, γd, γl, 0.03, 0.0, 0.0, 0.0, 0.0, d)) for d in Δ]; color = PAL[4]),
        _line("Hartmann-Tran", Δ,
              [real(pcqsdhc(cpf, 0.0, γd, γl, 0.03, 0.0, 0.0, 0.05, 0.4, d)) for d in Δ]; color = PAL[5]),
    ]
    layout = _layout("Area-normalized line-shape families (same γd, γl)", "ν − ν₀  [cm⁻¹]", "ϕ(ν)  [cm]")
    _write_plotly_asset("lineshape_families.html", data, layout; title = "Line-shape families")
end

# CO2 1.6 µm band at three temperatures — line strengths and widths shift with T.
function plot_temperature()
    par = joinpath(REFDATA, "golden", "co2_6300_6400_p500_T250.par")
    db = load_lines(HitranPort(par); mol = :CO2, iso = :main, min_strength = 1e-26)
    model = LineByLineModel(db; profile = Voigt(), wing_cutoff = 40.0)   # partition rides on db
    grid = collect(6300.0:0.02:6400.0)
    data = [_line("$(Int(T)) K", grid, clip(compute_cross_section(model, grid, 500.0, T)); color = c)
            for (T, c) in ((296.0, PAL[1]), (250.0, PAL[2]), (220.0, PAL[3]))]
    layout = _layout("CO₂ 1.6 µm band vs temperature (p = 500 hPa)", WAVENO, SIG; ylog = true)
    _write_plotly_asset("co2_temperature.html", data, layout; title = "CO₂ vs temperature")
end

# Read a committed cross-section overlay (cols: ν, σ_a, σ_b). Generated offline by the
# benchmark/gen_*.jl scripts so the doc build needs no network/HITRAN key (CI-safe).
function read_xsec(name)
    nu, a, b = Float64[], Float64[], Float64[]
    for ln in eachline(joinpath(REFDATA, "figures", name))
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
    data = [_line("ExoMol Li2015 (S from Einstein-A)", ν, clip(σe); color = PAL[1]),
            _line("HITRAN", ν, clip(σh); color = PAL[2], dash = "dot")]
    layout = _layout("CO cross-section — ExoMol vs HITRAN (p = 1013 hPa, T = 296 K)", WAVENO, SIG; ylog = true)
    _write_plotly_asset("exomol_co_xsec.html", data, layout; title = "ExoMol vs HITRAN CO")
end

# CO2 4.3 µm band: first-order line mixing redistributes intensity between overlapping
# lines and suppresses the far wings (sub-Lorentzian) vs a plain Voigt sum.
function plot_co2_linemix()
    ν, σno, σlm = read_xsec("co2_linemix_band.txt")
    data = [_line("Voigt (no line mixing)", ν, clip(σno); color = PAL[1]),
            _line("with line mixing (Hartmann-Tran)", ν, clip(σlm); color = PAL[2])]
    layout = _layout("CO₂ 4.3 µm band — effect of line mixing (p = 1013 hPa, T = 296 K)", WAVENO, SIG; ylog = true)
    _write_plotly_asset("co2_linemix.html", data, layout; title = "CO₂ line mixing")
end

# H2O: the speed dependence of collisions (Hartmann-Tran γ₂/ν_VC) narrows and reshapes the
# line cores relative to a plain Voigt.
function plot_h2o_ht()
    ν, σv, σht = read_xsec("h2o_voigt_vs_ht.txt")
    data = [_line("Voigt", ν, clip(σv); color = PAL[1]),
            _line("Hartmann-Tran (γ₂, ν_VC)", ν, clip(σht); color = PAL[2], dash = "dash")]
    layout = _layout("H₂O speed dependence — Voigt vs Hartmann-Tran (p = 1013 hPa, T = 250 K)", WAVENO, SIG; ylog = true)
    _write_plotly_asset("h2o_voigt_vs_ht.html", data, layout; title = "H₂O Voigt vs Hartmann-Tran")
end

function plot_benchmark()
    cut = ["2.5", "5", "10", "25", "full"]
    bench(name, y, c) = (; type = "scatter", mode = "lines+markers", name,
                         x = cut, y = collect(Float64, y), line = (; width = 2.5, color = c), marker = (; size = 8))
    data = [bench("hapi2 (numba)", [114.0, 173.0, 286.0, 623.0, 4691.0], PAL[2]),
            bench("AtmosphericAbsorption CPU", [177.0, 235.0, 349.0, 683.0, 4122.0], PAL[1]),
            bench("AtmosphericAbsorption GPU", [2.1, 2.2, 2.4, 3.0, 10.1], PAL[3])]
    layout = _layout("Time vs wing cutoff — 4000 lines × 400 cm⁻¹ band (A100)",
                     "wing cutoff  [cm⁻¹]", "time  [ms]"; ylog = true, xcategory = true)
    _write_plotly_asset("benchmark.html", data, layout; title = "Benchmark")
end

# Write every interactive figure into the built site. Data-free figures always render; the
# data-backed ones need the reference_data artifact. Called AFTER makedocs() (needs build/).
function generate_figure_assets()
    @info "Writing interactive figure assets into the built site…"
    plot_lineshape_families()      # pure math — no data needed
    plot_benchmark()               # tabulated timings — no data needed
    if REFDATA === nothing
        @warn "Skipping data-backed figures (temperature, CO₂ line-mixing, H₂O HT, ExoMol) — reference_data artifact unavailable"
    else
        plot_temperature()
        plot_co2_linemix()
        plot_h2o_ht()
        plot_exomol_co()
    end
end

# Generate the species/isotopologue reference table from the bundled data + registry, so
# the docs page stays in sync with the code that resolves :CO2 / :main / formulas. This is a
# SOURCE page, so it must run BEFORE makedocs().
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
        mm = fmt(() -> AtmosphericAbsorption.LineLists.molar_mass(r.mol, r.iso))
        ab = fmt(() -> AtmosphericAbsorption.LineLists.abundance(r.mol, r.iso))
        println(io, "| $(r.mol) | `:$sym` | $(r.iso) | `$(r.formula)` | $mm | $ab |")
    end
    write(joinpath(@__DIR__, "src", "isotopologues.md"), String(take!(io)))
end

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

# After the Vitepress build exists, write the interactive figures into it (and before deploy).
generate_figure_assets()

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
