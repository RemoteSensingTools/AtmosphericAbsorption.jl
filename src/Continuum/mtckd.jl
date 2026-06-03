#=
MT_CKD water-vapor continuum (self + foreign), MT_CKD v4.2 (Mlawer et al. 2012). The AER
reference coefficients ship bundled as data/mtckd.arrow (converted from the AER NetCDF).
The continuum cross-section follows the LBLRTM convention:

  σ_self(ν,T) = C_self(ν)·radterm(ν,T)·(p_h2o/p_ref)·(T_ref/T)^texp(ν)
  σ_for(ν,T)  = C_for(ν) ·radterm(ν,T)·(p_dry/p_ref)
  radterm(ν,T) = ν·tanh(c₂·ν/2T)

converting the AER coefficient [cm²/molec/cm⁻¹] into a cross-section [cm²/molecule]. The
per-layer τ = (σ_self+σ_for)·n_h2o·Δz is an RT-side concern. The table covers
ν ∈ [-20, 20000] cm⁻¹; outside, the continuum is zero.
=#

const _MTCKD_PATH = joinpath(@__DIR__, "..", "..", "data", "mtckd.arrow")
include_dependency(_MTCKD_PATH)

"""MT_CKD reference table: continuum coefficients at the native ν grid + reference state."""
struct MTCKDTable
    ν::Vector{Float64}            # cm⁻¹, ascending
    C_self::Vector{Float64}       # cm²/molec/cm⁻¹ at T_ref
    C_for::Vector{Float64}
    self_texp::Vector{Float64}
    p_ref::Float64                # hPa
    T_ref::Float64                # K
end

"""MT_CKD coefficients interpolated onto a model ν grid (zero outside the table range)."""
struct MTCKDBand
    ν::Vector{Float64}
    C_self::Vector{Float64}
    C_for::Vector{Float64}
    texp::Vector{Float64}
    p_ref::Float64
    T_ref::Float64
end

"""
    load_mtckd(path = bundled) -> MTCKDTable

Load MT_CKD reference coefficients from an Arrow file (defaults to the bundled AER v4.2
table). `p_ref`/`T_ref` come from the file's Arrow metadata.
"""
function load_mtckd(path::AbstractString = _MTCKD_PATH)
    t = Arrow.Table(path)
    md = Arrow.getmetadata(t)
    pref = md === nothing ? 1013.0 : parse(Float64, md["p_ref"])
    Tref = md === nothing ? 296.0 : parse(Float64, md["T_ref"])
    return MTCKDTable(collect(t.ν), collect(t.C_self), collect(t.C_for),
                      collect(t.self_texp), pref, Tref)
end

"""
    build_mtckd_band(table, ν_grid) -> MTCKDBand

Linearly interpolate the table onto `ν_grid`; coefficients are zero where `ν_grid` falls
outside the table's range (no continuum in the UV/Vis above ~500 nm).
"""
function build_mtckd_band(table::MTCKDTable, ν_grid::AbstractVector)
    n = length(ν_grid)
    Cs = zeros(n); Cf = zeros(n); te = zeros(n)
    νlo, νhi, N = table.ν[1], table.ν[end], length(table.ν)
    @inbounds for (k, νq_) in enumerate(ν_grid)
        νq = Float64(νq_)
        (νq < νlo || νq > νhi) && continue
        j = searchsortedfirst(table.ν, νq)
        if j ≤ 1
            Cs[k], Cf[k], te[k] = table.C_self[1], table.C_for[1], table.self_texp[1]
        elseif j > N
            Cs[k], Cf[k], te[k] = table.C_self[N], table.C_for[N], table.self_texp[N]
        else
            w = (νq - table.ν[j-1]) / (table.ν[j] - table.ν[j-1])
            Cs[k] = (1 - w) * table.C_self[j-1]    + w * table.C_self[j]
            Cf[k] = (1 - w) * table.C_for[j-1]     + w * table.C_for[j]
            te[k] = (1 - w) * table.self_texp[j-1] + w * table.self_texp[j]
        end
    end
    return MTCKDBand(Float64.(ν_grid), Cs, Cf, te, table.p_ref, table.T_ref)
end

"""
    h2o_continuum!(σ_out, band, T, p_h2o, p_dry) -> σ_out

Total (self + foreign) H₂O continuum cross-section [cm²/molecule] on `band.ν` at
temperature `T` [K] with H₂O and dry-air partial pressures `p_h2o`, `p_dry` [hPa].
"""
function h2o_continuum!(σ_out::AbstractVector{Float64}, band::MTCKDBand,
                        T::Real, p_h2o::Real, p_dry::Real)
    c₂ = C2_RAD
    self_p, for_p = p_h2o / band.p_ref, p_dry / band.p_ref
    @inbounds for k in eachindex(band.ν)
        ν = band.ν[k]
        radterm = ν * tanh(c₂ * ν / (2 * T))
        σ_out[k] = band.C_self[k] * radterm * self_p * (band.T_ref / T)^band.texp[k] +
                   band.C_for[k]  * radterm * for_p
    end
    return σ_out
end

"""
    h2o_continuum(band, T, p_h2o, p_dry) -> Vector{Float64}

Allocating form of [`h2o_continuum!`](@ref).
"""
h2o_continuum(band::MTCKDBand, T::Real, p_h2o::Real, p_dry::Real) =
    h2o_continuum!(Vector{Float64}(undef, length(band.ν)), band, T, p_h2o, p_dry)
