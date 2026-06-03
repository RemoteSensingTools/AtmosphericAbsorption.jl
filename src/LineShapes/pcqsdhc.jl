#=
The pCqSDHC line shape (partially-Correlated quadratic-Speed-Dependent Hard-Collision,
Tran/Ngo/Hartmann, JQSRT 2013) — the engine behind the Hartmann-Tran (HTP) and
speed-dependent Voigt profiles. Ported from HAPI's `pcqsdhc`. Returns the FULL complex
normalized shape; the real part is the profile, the imaginary part feeds line mixing.

Parameters (all cm⁻¹ except η): Γ0/Δ0 speed-averaged width/shift, Γ2/Δ2 their speed
dependence, νVC velocity-changing frequency, η correlation. With Γ2=Δ2=νVC=η=0 this
reduces exactly to the Voigt profile.
=#

"""
    pcqsdhc(cpf, ν0, γd, Γ0, Γ2, Δ0, Δ2, νVC, η, ν) -> Complex

Complex normalized pCqSDHC line shape at wavenumber `ν` for a line at `ν0` with Doppler
HWHM `γd`. `real` is the line profile (cm), `imag` is used by first-order line mixing.
"""
# Shared A/B-term pattern of PART 2 and PART 4 (two CPF evaluations at Z1, Z2).
@inline _aterm(rpi, cte, W1, W2) = rpi * cte * (W1 - W2)
@inline function _bterm(rpi, csqrtY, c2t, Z1, W1, Z2, W2)
    h = rpi / (2 * csqrtY)
    return (-1 + h * (1 - Z1^2) * W1 - h * (1 - Z2^2) * W2) / c2t
end

@inline function pcqsdhc(cpf::AbstractCPF, ν0::FT, γd::FT, Γ0::FT, Γ2::FT,
                         Δ0::FT, Δ2::FT, νVC::FT, η::FT, ν::FT) where {FT}
    cte = FT(SQRT_LN2) / γd
    rpi = FT(SQRT_PI)
    c0  = complex(Γ0, Δ0)
    c2  = complex(Γ2, Δ2)
    c0t = (1 - η) * (c0 - FT(1.5) * c2) + νVC
    c2t = (1 - η) * c2
    Δ   = ν0 - ν

    if iszero(c2t)                                        # PART 1: no speed-dependent width
        Z1 = (im * Δ + c0t) * cte
        W1 = w(cpf, im * Z1)
        Aterm = rpi * cte * W1
        Bterm = abs(Z1) ≤ 4000 ?
            rpi * cte * ((1 - Z1^2) * W1 + Z1 / rpi) :
            cte * (rpi * W1 + FT(0.5) / Z1 - FT(0.75) / Z1^3)
    else
        X      = (im * Δ + c0t) / c2t
        Y      = inv((2 * cte * c2t)^2)
        csqrtY = (Γ2 - im * Δ2) / (2 * cte * (1 - η) * (Γ2^2 + Δ2^2))
        if abs(X) ≤ FT(3e-8) * abs(Y)                     # PART 2
            Z1 = (im * Δ + c0t) * cte
            Z2 = sqrt(X + Y) + csqrtY
            W1, W2 = w(cpf, im * Z1), w(cpf, im * Z2)
            Aterm = _aterm(rpi, cte, W1, W2)
            Bterm = _bterm(rpi, csqrtY, c2t, Z1, W1, Z2, W2)
        elseif abs(Y) ≤ FT(1e-15) * abs(X)                # PART 3: large-c2t asymptotic (unphysical for real spectra)
            sqrtXY = sqrt(X + Y)
            W1 = w(cpf, im * sqrtXY)
            if abs(sqrt(X)) ≤ 4000
                sqrtX = sqrt(X)
                Wb = w(cpf, im * sqrtX)
                Aterm = (2rpi / c2t) * (1 / rpi - sqrtX * Wb)
                Bterm = (1 / c2t) * (-1 + 2rpi * (1 - X - 2Y) * (1 / rpi - sqrtX * Wb)
                                        + 2rpi * sqrtXY * W1)
            else
                Aterm = (1 / c2t) * (1 / X - FT(1.5) / X^2)
                Bterm = (1 / c2t) * (-1 + (1 - X - 2Y) * (1 / X - FT(1.5) / X^2)
                                        + 2rpi * sqrtXY * W1)
            end
        else                                              # PART 4 (the general case)
            Z1 = sqrt(X + Y) - csqrtY
            Z2 = Z1 + 2 * csqrtY
            W1, W2 = w(cpf, im * Z1), w(cpf, im * Z2)
            Aterm = _aterm(rpi, cte, W1, W2)
            Bterm = _bterm(rpi, csqrtY, c2t, Z1, W1, Z2, W2)
        end
    end

    return Aterm / (FT(π) * (1 - (νVC - η * (c0 - FT(1.5) * c2)) * Aterm + η * c2 * Bterm))
end

# Uniform per-line entry points (see profiles.jl). SDV is pCqSDHC with νVC=η=0. First-
# order (Rosenkranz) line mixing folds in as real(LS) + Y·imag(LS), reusing the complex
# shape — free, and zero when the line list carries no mixing coefficient (Y=0).
@inline function evaluate(::SpeedDependentVoigt, cpf::AbstractCPF, νI::FT, ν0::FT, p) where {FT}
    LS = pcqsdhc(cpf, ν0, p.γd, p.Γ0, p.Γ2, p.Δ0, p.Δ2, zero(FT), zero(FT), νI)
    return real(LS) + p.Y * imag(LS)
end
@inline function evaluate(::HartmannTran, cpf::AbstractCPF, νI::FT, ν0::FT, p) where {FT}
    LS = pcqsdhc(cpf, ν0, p.γd, p.Γ0, p.Γ2, p.Δ0, p.Δ2, p.νVC, p.η, νI)
    return real(LS) + p.Y * imag(LS)
end
