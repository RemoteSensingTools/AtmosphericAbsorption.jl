"""
    SourceMetadata

Provenance and reference state for a line list, carried alongside the columns.
"""
struct SourceMetadata
    source::String     # e.g. "HITRAN2024 CO2"
    T_ref::Float64     # reference temperature [K]
    p_ref::Float64     # reference pressure [hPa]
end

"""
    LineDatabase{FT}

Struct-of-arrays line list holding the union of parameters every profile needs.
Core Voigt columns are always present; advanced columns (speed-dependence, velocity
changing, line mixing) default to zero so one layout and one kernel serve every
profile. All float columns share element type `FT`; all share length `length(db)`.

The list is host-resident (built on CPU by a Port); only the prepared per-line
parameters are uploaded to the device. `ν0` is expected sorted ascending so Ports
and the pre-pass can window with binary search.

Construct with the keyword builder `LineDatabase(; mol, iso, ν0, …)` — advanced
columns omitted are zero-filled.
"""
struct LineDatabase{FT<:AbstractFloat}
    mol::Vector{Int32}      # HITRAN molecule id
    iso::Vector{Int32}      # isotopologue id
    ν0::Vector{FT}          # line center [cm⁻¹]
    S::Vector{FT}           # intensity at T_ref [cm⁻¹/(molecule·cm⁻²)]
    E_lower::Vector{FT}     # lower-state energy [cm⁻¹]
    g_upper::Vector{FT}     # upper-state degeneracy
    γ_air::Vector{FT}       # air-broadened HWHM [cm⁻¹/atm] at T_ref
    γ_self::Vector{FT}      # self-broadened HWHM [cm⁻¹/atm]
    n_air::Vector{FT}       # T-exponent of γ_air
    δ_air::Vector{FT}       # air pressure shift [cm⁻¹/atm]
    n_self::Vector{FT}      # T-exponent of γ_self (defaults to n_air)
    δ_self::Vector{FT}      # self pressure shift [cm⁻¹/atm] (defaults to 0)
    γ2_air::Vector{FT}      # speed-dependent broadening (HT/SDV)
    δ2_air::Vector{FT}      # speed-dependent shift     (HT/SDV)
    νVC::Vector{FT}         # velocity-changing collision frequency (HT)
    η::Vector{FT}           # correlation parameter (HT)
    n_γ2::Vector{FT}        # T-exponent of γ2
    Y_LM::Vector{FT}        # first-order line-mixing coefficient at T_ref
    n_Y_LM::Vector{FT}      # T-exponent of Y_LM
    molar_mass::Vector{FT}  # per-line molar mass [g/mol] (for the Doppler width)
    meta::SourceMetadata
end

Base.length(db::LineDatabase) = length(db.ν0)
Base.eltype(::LineDatabase{FT}) where {FT} = FT

function LineDatabase(; mol, iso, ν0, S, E_lower, g_upper,
                        γ_air, γ_self, n_air, δ_air, molar_mass, meta,
                        n_self = nothing, δ_self = nothing,
                        γ2_air = nothing, δ2_air = nothing, νVC = nothing,
                        η = nothing, n_γ2 = nothing, Y_LM = nothing, n_Y_LM = nothing)
    FT = eltype(ν0)
    n  = length(ν0)
    col(x) = x === nothing ? zeros(FT, n) : convert(Vector{FT}, x)
    nair = convert(Vector{FT}, n_air)
    # No self exponent ⇒ scale the self width with the air exponent (HAPI's fallback);
    # no self shift ⇒ zero.
    nself = n_self === nothing ? copy(nair) : convert(Vector{FT}, n_self)
    return LineDatabase(convert(Vector{Int32}, mol), convert(Vector{Int32}, iso),
                        ν0, S, E_lower, g_upper, γ_air, γ_self, nair, δ_air,
                        nself, col(δ_self),
                        col(γ2_air), col(δ2_air), col(νVC), col(η), col(n_γ2),
                        col(Y_LM), col(n_Y_LM), convert(Vector{FT}, molar_mass), meta)
end
