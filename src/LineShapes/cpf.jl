#=
Complex probability function w(z) = exp(-z²)·erfc(-iz), the kernel of the Voigt
profile. We expose ONE production strategy (`HumlicekWeideman32`, GPU-proven) plus
`ErfcxCPF` as a CPU-only reference for validation. `w` returns the FULL complex value
so line-mixing (which needs Im w) can reuse the same evaluation.
=#

"""Strategy for evaluating the complex probability function `w(z)`. Internal detail."""
abstract type AbstractCPF end

"""
    HumlicekWeideman32()

Region-split `w(z)`: Humlíček (1982) region-II rational for |x|+y ≥ 8, Weideman (1994)
32-term rational otherwise. Allocation-free and GPU-safe — the default everywhere.
"""
struct HumlicekWeideman32 <: AbstractCPF end

"""
    ErfcxCPF()

Reference `w(z) = erfcx(-iz)` via SpecialFunctions. Accurate but **CPU-only**
(scalar special function) — used to validate `HumlicekWeideman32`.
"""
struct ErfcxCPF <: AbstractCPF end

@inline w(::ErfcxCPF, z::Complex{FT}) where {FT} = erfcx(-im * z)

@inline function w(::HumlicekWeideman32, z::Complex{FT}) where {FT}
    return abs(real(z)) + imag(z) ≥ 8 ? _humlicek_region2(z) : _weideman32(z)
end

"""Humlíček (1982) region-II rational approximation (large |x|+y)."""
@inline function _humlicek_region2(z::Complex{FT}) where {FT}
    t = imag(z) - im * real(z)
    u = t^2
    return (t * (FT(1.410474) + u * FT(1 / sqrt(π)))) / (FT(0.75) + u * (3 + u))
end

"""Weideman (1994) 32-term rational approximation, SIAM J. Numer. Anal. 31, 1497."""
@inline function _weideman32(z::Complex{FT}) where {FT}
    L = FT(sqrt(32 / sqrt(2)))
    a = (FT(2.5722534081245696), FT(2.2635372999002676), FT(1.8256696296324824),
         FT(1.3455441692345453), FT(9.0192548936480144e-1), FT(5.4601397206393498e-1),
         FT(2.9544451071508926e-1), FT(1.4060716226893769e-1), FT(5.7304403529837900e-2),
         FT(1.9006155784845689e-2), FT(4.5195411053501429e-3), FT(3.9259136070122748e-4),
         FT(-2.4532980269928922e-4), FT(-1.3075449254548613e-4), FT(-2.1409619200870880e-5),
         FT(6.8210319440412389e-6), FT(4.4015317319048931e-6), FT(4.2558331390536872e-7),
         FT(-4.1840763666294341e-7), FT(-1.4813078891201116e-7), FT(2.2930439569075392e-8),
         FT(2.3797557105844622e-8), FT(8.1248960947953431e-10), FT(-3.2080150458594088e-9),
         FT(-5.2310170266050247e-10), FT(4.1537465934749353e-10), FT(1.1658312885903929e-10),
         FT(-5.5441820344468828e-11), FT(-2.1542618451370239e-11), FT(8.0314997274316680e-12),
         FT(3.7424975634801558e-12), FT(-1.3031797863050087e-12))

    iz       = im * real(z) - imag(z)
    rec_lmiz = inv(L - iz)
    Z        = (L + iz) * rec_lmiz
    p = complex(a[32])                  # Complex{FT} from the start: keeps the Horner loop GPU-safe
    for k in 31:-1:1
        @inbounds p = a[k] + p * Z
    end
    return (FT(1 / sqrt(π)) + 2 * p * rec_lmiz) * rec_lmiz
end
