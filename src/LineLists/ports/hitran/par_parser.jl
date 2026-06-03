#=
Parser for the HITRAN 160-character fixed-width line-by-line (`.par`) format. We read
only the fields the cross-section needs, filtering by molecule/isotopologue/window/
strength during the scan so huge files never fully materialize.
=#

# 1-indexed character ranges of the HITRAN .par fields we use.
const _PAR = (mol = 1:2, iso = 3:3, ν = 4:15, S = 16:25, γ_air = 36:40,
              γ_self = 41:45, E = 46:55, n_air = 56:59, δ = 60:67, g_up = 147:153)

@inline _f64(line, r) = parse(Float64, @view line[r])

# HITRAN local_iso_id: digits 1–9, '0' → 10, then 'A' → 11, 'B' → 12, …
@inline function _isoid(c::Char)
    c == '0' && return 10
    isdigit(c) && return Int(c - '0')
    return Int(c - 'A') + 11
end

"Raw, column-oriented HITRAN line parameters straight from a `.par` file."
struct ParColumns
    mol::Vector{Int32}; iso::Vector{Int32}
    ν::Vector{Float64}; S::Vector{Float64}; E::Vector{Float64}
    γ_air::Vector{Float64}; γ_self::Vector{Float64}
    n_air::Vector{Float64}; δ::Vector{Float64}; g_up::Vector{Float64}
end

"""
    parse_par(path; mol=-1, iso=-1, ν_min=0.0, ν_max=Inf, min_strength=0.0) -> ParColumns

Scan a HITRAN `.par` file, keeping lines that match `mol`/`iso` (`-1` = any), fall in
`[ν_min, ν_max]`, and have intensity ≥ `min_strength`. Returns the raw columns in
file order.
"""
function parse_par(path::AbstractString; mol::Integer = -1, iso::Integer = -1,
                   ν_min::Real = 0.0, ν_max::Real = Inf, min_strength::Real = 0.0)
    cols = ParColumns(Int32[], Int32[], Float64[], Float64[], Float64[],
                      Float64[], Float64[], Float64[], Float64[], Float64[])
    for line in eachline(path)
        length(line) < 160 && continue
        m = parse(Int32, @view line[_PAR.mol])
        (mol == -1 || m == mol) || continue
        i = Int32(_isoid(line[3]))
        (iso == -1 || i == iso) || continue
        ν = _f64(line, _PAR.ν)
        ν_min ≤ ν ≤ ν_max || continue
        S = _f64(line, _PAR.S)
        S ≥ min_strength || continue
        push!(cols.mol, m);            push!(cols.iso, i)
        push!(cols.ν, ν);              push!(cols.S, S)
        push!(cols.E, _f64(line, _PAR.E))
        push!(cols.γ_air, _f64(line, _PAR.γ_air))
        push!(cols.γ_self, _f64(line, _PAR.γ_self))
        push!(cols.n_air, _f64(line, _PAR.n_air))
        push!(cols.δ, _f64(line, _PAR.δ))
        push!(cols.g_up, _f64(line, _PAR.g_up))
    end
    return cols
end
