#=
CPU pre-pass: turn a LineDatabase + (p, T) into the per-active-line scalars the
kernel needs (shifted center, Doppler/Lorentz widths, width ratio, T-corrected
intensity, grid index window), then upload them to the compute device. All work is
in `FT`; nothing promotes to Float64 on a Float32 model.
=#

"Prepared, device-resident line parameters for one (grid, p, T) evaluation."
struct PreparedLines{V,VI}
    ν::V; γd::V; γl::V; y::V; S::V    # shifted center, widths, ratio, intensity
    istart::VI; istop::VI            # inclusive grid index window per line
    n::Int
end

function prepare(model::LineByLineModel{FT}, grid::AbstractVector,
                 pressure::Real, temperature::Real) where {FT}
    (; lines, partition, wing_cutoff, vmr, architecture) = model
    T, p   = FT(temperature), FT(pressure)
    Tref   = FT(lines.meta.T_ref)
    pratio = p / FT(lines.meta.p_ref)
    c₂     = FT(C2_RAD)
    # ν₀-independent Doppler factor: γd = ν₀ · dopp / √(molar_mass[g/mol]),
    # since a molecule of mass M g/mol weighs M·AMU kg.
    dopp   = FT(SQRT_2LN2) / FT(C_LIGHT) * sqrt(FT(K_BOLTZMANN) * T / FT(AMU))

    # Select on the shifted center νs = ν₀ + pratio·δ_air. We pad the ν₀ band by the
    # largest possible shift so no edge line that shifts into the window is dropped;
    # the exact per-line window (istart/istop) is computed from νs below.
    δmax = isempty(lines.δ_air) ? zero(FT) : abs(pratio) * maximum(abs, lines.δ_air)
    νmin = FT(first(grid)) - wing_cutoff - δmax
    νmax = FT(last(grid)) + wing_cutoff + δmax
    active = findall(ν0 -> νmin ≤ ν0 ≤ νmax, lines.ν0)
    n = length(active)

    ν  = Vector{FT}(undef, n); γd = similar(ν); γl = similar(ν)
    y  = similar(ν);           S  = similar(ν)
    istart = Vector{Int32}(undef, n); istop = similar(istart)
    Qcache = Dict{Tuple{Int,Int},FT}()
    Ng = length(grid)

    @inbounds for (k, j) in enumerate(active)
        ν0 = lines.ν0[j]
        νs = ν0 + pratio * lines.δ_air[j]                                  # pressure-shifted center
        γl_j = (lines.γ_air[j] * (1 - vmr) + lines.γ_self[j] * vmr) *
               pratio * (Tref / T)^lines.n_air[j]
        γd_j = ν0 * dopp / sqrt(lines.molar_mass[j])
        key  = (Int(lines.mol[j]), Int(lines.iso[j]))
        Qr   = get(Qcache, key, zero(FT))
        if Qr == 0
            Qr = FT(Q_ratio(partition, lines.mol[j], lines.iso[j], T, Tref))
            Qcache[key] = Qr
        end
        S_j = lines.S[j] * Qr * exp(c₂ * lines.E_lower[j] * (1 / Tref - 1 / T)) *
              (-expm1(-c₂ * ν0 / T)) / (-expm1(-c₂ * ν0 / Tref))           # intensity at T

        ν[k] = νs; γl[k] = γl_j; γd[k] = γd_j
        y[k] = FT(SQRT_LN2) * γl_j / γd_j; S[k] = S_j
        istart[k] = clamp(searchsortedfirst(grid, νs - wing_cutoff), 1, Ng)
        istop[k]  = clamp(searchsortedlast(grid, νs + wing_cutoff), 1, Ng)
    end

    to = array_type(architecture)
    return PreparedLines(to(ν), to(γd), to(γl), to(y), to(S), to(istart), to(istop), n)
end
