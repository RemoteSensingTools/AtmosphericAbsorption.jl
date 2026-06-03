#=
CPU pre-pass: turn a LineDatabase + (p, T) into the per-active-line parameters every
profile needs — line center, Doppler width, the (speed-averaged and speed-dependent)
widths/shifts, the velocity-changing frequency and correlation, the T-corrected
intensity, and the grid index window — then upload them to the compute device. All
work is in `FT`; nothing promotes to Float64 on a Float32 model.
=#

"Prepared, device-resident line parameters for one (grid, p, T) evaluation."
struct PreparedLines{V,VI}
    ν0::V; γd::V                     # line center, Doppler HWHM
    Γ0::V; Γ2::V; Δ0::V; Δ2::V       # width, speed-dep width, shift, speed-dep shift
    νVC::V; η::V                     # velocity-changing freq, correlation (HT)
    S::V                            # T-corrected intensity
    istart::VI; istop::VI           # inclusive grid index window per line
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

    # Select on the shifted center νs = ν₀ + Δ0. We pad the ν₀ band by the largest
    # possible shift so no edge line that shifts into the window is dropped; the exact
    # per-line window (istart/istop) is computed from νs below. `ν0` is sorted
    # (LineDatabase contract), so the band is found with two binary searches.
    δmax = isempty(lines.δ_air) ? zero(FT) : abs(pratio) * maximum(abs, lines.δ_air)
    νmin = FT(first(grid)) - wing_cutoff - δmax
    νmax = FT(last(grid)) + wing_cutoff + δmax
    active = searchsortedfirst(lines.ν0, νmin):searchsortedlast(lines.ν0, νmax)
    n = length(active)

    ν0 = Vector{FT}(undef, n); γd = similar(ν0)
    Γ0 = similar(ν0); Γ2 = similar(ν0); Δ0 = similar(ν0); Δ2 = similar(ν0)
    νVC = similar(ν0); η = similar(ν0); S = similar(ν0)
    istart = Vector{Int32}(undef, n); istop = similar(istart)
    Qcache = Dict{Tuple{Int32,Int32},FT}()   # partition ratio per (mol, iso)
    Ng = length(grid)

    @inbounds for (k, j) in enumerate(active)
        ν0j = lines.ν0[j]
        Δ0j = pratio * lines.δ_air[j]                                       # pressure shift
        Γ0j = (lines.γ_air[j] * (1 - vmr) + lines.γ_self[j] * vmr) *
              pratio * (Tref / T)^lines.n_air[j]                           # speed-averaged width
        key = (lines.mol[j], lines.iso[j])
        Qr  = get(Qcache, key, FT(NaN))             # NaN sentinel: never a valid ratio
        if isnan(Qr)
            Qr = FT(Q_ratio(partition, key[1], key[2], T, Tref))
            Qcache[key] = Qr
        end
        S[k]  = lines.S[j] * Qr * exp(c₂ * lines.E_lower[j] * (1 / Tref - 1 / T)) *
                (-expm1(-c₂ * ν0j / T)) / (-expm1(-c₂ * ν0j / Tref))        # intensity at T
        ν0[k] = ν0j
        γd[k] = ν0j * dopp / sqrt(lines.molar_mass[j])
        Γ0[k] = Γ0j;                              Δ0[k] = Δ0j
        Γ2[k] = lines.γ2_air[j] * pratio * (Tref / T)^lines.n_γ2[j]        # advanced cols are
        Δ2[k] = lines.δ2_air[j] * pratio                                   # zero for Voigt lists
        νVC[k] = lines.νVC[j] * pratio;          η[k] = lines.η[j]
        νs = ν0j + Δ0j
        istart[k] = clamp(searchsortedfirst(grid, νs - wing_cutoff), 1, Ng)
        istop[k]  = clamp(searchsortedlast(grid, νs + wing_cutoff), 1, Ng)
    end

    to = array_type(architecture)
    return PreparedLines(to(ν0), to(γd), to(Γ0), to(Γ2), to(Δ0), to(Δ2),
                         to(νVC), to(η), to(S), to(istart), to(istop), n)
end
