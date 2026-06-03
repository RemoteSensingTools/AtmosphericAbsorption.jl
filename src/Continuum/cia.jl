#=
HITRAN Collision-Induced Absorption (CIA). Reads HITRAN `.cia` files ((T, ν, σ) blocks
per pair, e.g. O2-O2, O2-N2, N2-N2), projects σ onto a model ν grid per block
temperature, and interpolates in temperature. Exposes the binary cross-section
σ_AB(ν, T) [cm⁵/molecule²]; the per-layer optical depth τ = σ_AB·n_A·n_B·Δz is an
RT-side concern. All arithmetic is Float64 (σ ~1e-44..1e-46, below Float32's range).
=#

"One (T, ν, σ) block from a HITRAN `.cia` file."
struct CIABlock{FT}
    formula::String
    T::FT
    ν::Vector{FT}
    σ::Vector{FT}     # cm⁵/molecule²
end

"""
    CIATable

CIA cross-section with each block projected onto a ν grid. `σ_grid` is
`(length(ν_grid), n_blocks)`, holding block `b`'s σ where it covers the grid point and
`NaN` elsewhere; `block_T` are the (ascending) block temperatures. Temperature
interpolation is done **per frequency** over only the blocks covering it
([`cia_cross_section!`](@ref)), so windows spanning several bands measured at different
temperatures stay correct.
"""
struct CIATable
    pair::String                        # e.g. "O2-O2"
    species_a::String
    species_b::String
    ν::Vector{Float64}                  # model grid
    block_T::Vector{Float64}            # per-block temperature, ascending
    σ_grid::Matrix{Float64}             # (length(ν), n_blocks); NaN outside a block's ν range
end

"""
    parse_cia_file(path) -> Vector{CIABlock{Float64}}

Read a HITRAN `.cia` file: each block is a header line then `n_pts` `(ν σ)` rows. Header
(fixed columns): formula [1:20], n_pts [41:47], T [48:54].
"""
function parse_cia_file(path::AbstractString)
    blocks = CIABlock{Float64}[]
    open(path) do io
        while !eof(io)
            header = readline(io)
            (length(header) < 54 || isempty(strip(header))) && continue
            formula = String(strip(header[1:20]))
            n_pts   = parse(Int, strip(header[41:47]))
            T_K     = parse(Float64, strip(header[48:54]))
            νs = Vector{Float64}(undef, n_pts); σs = Vector{Float64}(undef, n_pts)
            for i in 1:n_pts
                vals = split(strip(readline(io)))
                νs[i] = parse(Float64, vals[1]); σs[i] = parse(Float64, vals[2])
            end
            push!(blocks, CIABlock{Float64}(formula, T_K, νs, σs))
        end
    end
    return blocks
end

# Linear-interpolate one block onto ν_grid (in-place); leave NaN outside the block's range.
function _project_block!(σ_out, ν_grid, ν_blk, σ_blk)
    νlo, νhi, n = ν_blk[1], ν_blk[end], length(ν_blk)
    for (k, νq) in enumerate(ν_grid)
        (νq < νlo || νq > νhi) && continue
        j = searchsortedfirst(ν_blk, νq)
        σ_out[k] = j ≤ 1 ? σ_blk[1] : j > n ? σ_blk[end] :
                   let w = (νq - ν_blk[j-1]) / (ν_blk[j] - ν_blk[j-1])
                       (1 - w) * σ_blk[j-1] + w * σ_blk[j]
                   end
    end
    return σ_out
end

_split_pair(formula) = (p = split(strip(formula), '-'); length(p) ≥ 2 ?
    (String(p[1]), String(p[2])) : error("CIA pair \"$formula\" not \"A-B\""))

"""
    build_cia_table(blocks, ν_grid) -> CIATable

Project each window-overlapping block onto `ν_grid` (NaN where a block doesn't cover a
grid point), keeping blocks sorted by temperature so temperature interpolation can be done
per frequency over only the covering blocks.
"""
function build_cia_table(blocks::AbstractVector{<:CIABlock}, ν_grid::AbstractVector)
    isempty(blocks) && error("build_cia_table: no blocks")
    νf = Float64.(ν_grid)
    gmin, gmax = extrema(νf)
    relevant = filter(blk -> blk.ν[end] ≥ gmin && blk.ν[1] ≤ gmax, blocks)
    a, b = _split_pair(blocks[1].formula)
    perm = sortperm(Float64[blk.T for blk in relevant])
    relevant = relevant[perm]
    σg = fill(NaN, length(νf), length(relevant))
    for (jb, blk) in enumerate(relevant)
        _project_block!(view(σg, :, jb), νf, blk.ν, blk.σ)
    end
    return CIATable(blocks[1].formula, a, b, νf, Float64[blk.T for blk in relevant], σg)
end

"""
    load_cia(path, ν_grid) -> CIATable

Parse a `.cia` file and project it onto `ν_grid`.
"""
load_cia(path::AbstractString, ν_grid::AbstractVector) = build_cia_table(parse_cia_file(path), ν_grid)

# Interpolate one frequency's σ in temperature over the blocks covering it (NaN = not
# covered). Linear between bracketing covering blocks; constant extrapolation; 0 if none.
@inline function _interp_T(σrow, Ts, T)
    lo = 0
    @inbounds for b in eachindex(Ts)
        isnan(σrow[b]) && continue
        if Ts[b] ≤ T
            lo = b
        else
            lo == 0 && return σrow[b]                         # below all → lowest
            w = (T - Ts[lo]) / (Ts[b] - Ts[lo])
            return (1 - w) * σrow[lo] + w * σrow[b]
        end
    end
    return lo == 0 ? 0.0 : σrow[lo]                           # above all → highest, or none → 0
end

"""
    cia_cross_section!(σ_out, table, T) -> σ_out

Binary CIA cross-section σ(ν, T) [cm⁵/molecule²] on the table's ν grid. Temperature
interpolation is per frequency over only the blocks covering it (linear; constant
extrapolation; zero where no block covers the frequency).
"""
function cia_cross_section!(σ_out::AbstractVector{Float64}, table::CIATable, T::Real)
    Tf, Ts = Float64(T), table.block_T
    @inbounds for k in eachindex(σ_out)
        σ_out[k] = isempty(Ts) ? 0.0 : _interp_T(view(table.σ_grid, k, :), Ts, Tf)
    end
    return σ_out
end

"""
    cia_cross_section(table, T) -> Vector{Float64}

Allocating form of [`cia_cross_section!`](@ref).
"""
cia_cross_section(table::CIATable, T::Real) =
    cia_cross_section!(Vector{Float64}(undef, length(table.ν)), table, T)
