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
columns omitted are zero-filled. The partition function `Q(T)` for the source rides
along (`partition`; set by the Port), so a model needs no separate one.
"""
struct LineDatabase{FT<:AbstractFloat,PF<:AbstractPartitionFunction}
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
    partition::PF           # partition function Q(T) for the source (TIPS, ExoMol .pf)
end

Base.length(db::LineDatabase) = length(db.ν0)
Base.eltype(::LineDatabase{FT}) where {FT} = FT

function LineDatabase(; mol, iso, ν0, S, E_lower, g_upper,
                        γ_air, γ_self, n_air, δ_air, molar_mass, meta,
                        partition = TIPS2021PF(),
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
                        col(Y_LM), col(n_Y_LM), convert(Vector{FT}, molar_mass), meta, partition)
end

"""
    molecules(db) -> Vector{Symbol}

The distinct molecules present in `db`, as formula symbols (e.g. `[:H2O, :CO2]`) — the
"which species did this band pick up?" query.
"""
molecules(db::LineDatabase) = [molecule_symbol(m) for m in sort(unique(db.mol))]

"""
    db[mask]   -> LineDatabase
    db[:CO2]   -> LineDatabase

Subset a line list by a boolean mask (`db[db.mol .== 2]`) or by molecule (`db[:CO2]` /
`db["CO2"]`), preserving the partition function and provenance — e.g. to split a
multi-molecule band and plot each species independently.
"""
function Base.getindex(db::LineDatabase, mask::AbstractVector{Bool})
    LineDatabase(; mol = db.mol[mask], iso = db.iso[mask], ν0 = db.ν0[mask], S = db.S[mask],
        E_lower = db.E_lower[mask], g_upper = db.g_upper[mask], γ_air = db.γ_air[mask],
        γ_self = db.γ_self[mask], n_air = db.n_air[mask], δ_air = db.δ_air[mask],
        n_self = db.n_self[mask], δ_self = db.δ_self[mask], γ2_air = db.γ2_air[mask],
        δ2_air = db.δ2_air[mask], νVC = db.νVC[mask], η = db.η[mask], n_γ2 = db.n_γ2[mask],
        Y_LM = db.Y_LM[mask], n_Y_LM = db.n_Y_LM[mask], molar_mass = db.molar_mass[mask],
        meta = db.meta, partition = db.partition)
end
Base.getindex(db::LineDatabase, m::Union{Symbol,AbstractString}) =
    db[db.mol .== molecule_number(m)]

# Distinct (mol, iso) pairs, as "CO2 (12C)(16O)2" labels, for display.
function _species_labels(db::LineDatabase)
    pairs = sort!(collect(Set(zip(db.mol, db.iso))))   # dedup without materializing the full zip
    return ["$(molecule_symbol(m)) $(isotopologue(m, i))" for (m, i) in pairs]
end

function Base.show(io::IO, db::LineDatabase{FT}) where {FT}
    print(io, "LineDatabase{", FT, "}(", length(db), " lines")
    isempty(db.mol) || print(io, ", ", join(molecules(db), ", "))
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", db::LineDatabase{FT}) where {FT}
    n = length(db)
    print(io, "LineDatabase{", FT, "} — ", n, n == 1 ? " transition" : " transitions")
    if n > 0
        labels = _species_labels(db)
        print(io, "\n  species:    ", length(labels) ≤ 6 ? join(labels, ", ") :
              "$(length(labels)) isotopologues across $(join(molecules(db), ", "))")
        @printf(io, "\n  range:      %.5g – %.5g cm⁻¹", first(db.ν0), last(db.ν0))
    end
    print(io, "\n  source:     ", db.meta.source)
    print(io, "\n  partition:  ", pf_name(db.partition))
end
