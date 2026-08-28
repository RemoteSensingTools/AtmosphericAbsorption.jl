# Batched cross-section queries: many (p, T[, broadener]) points in ONE kernel.
#
# WHY. The tabulated models (InterpolationModel, AbscoLUT) are architecture-
# aware, but their scalar query is a poor fit for a GPU when a caller needs a
# whole atmospheric profile: vSmartMOM's `compute_absorption_profile!` used to
# loop 72 layers on the host, and each layer's query cost ~7 broadcast kernel
# launches (T/vmr/p blends, each materializing a table-window temporary), a
# host->device upload of the query grid, a FULL device synchronize, and a
# device->host copy of the result. Per model build that is ~700 launches, 144
# transfers and 72 pipeline drains for microseconds of arithmetic each — the
# GPU sits idle while launch latency dominates.
#
# The batched entry point takes VECTORS of query points and evaluates the
# entire (grid x point) surface in one kernel: brackets are precomputed on the
# host (they live on host axes anyway), uploaded once, and each work item
# fuses the (p, T[, vmr]) blend with the nu-resample for one output element.
# One launch, one synchronize, and the result STAYS on the architecture.
#
# The point axis is deliberately anonymous: a caller may pass 72 layers of one
# atmospheric column, or the concatenated layers of MANY columns (multi-cell
# batching) — the kernel only sees a flat list of (p, T, vmr) queries, so
# efficiency scales with total points, not with how the caller organizes them.
#
# NUMERICS. Each work item evaluates the SAME expression tree, in the same
# order, as the scalar path it replaces: T-blend, then vmr-blend, then
# p-blend, then the nu lerp, with identical `_bracket` semantics and identical
# zero-outside-range behavior. On the CPU backend the results are bit-identical
# to the scalar path; a GPU may fuse multiply-add differently, which is the
# same caveat every other kernel in this package carries.

"""
    compute_cross_section_profile(model, grid, pressures, temperatures;
                                  vmr = nothing, interp = :linear) -> matrix

Cross-sections [cm²/molecule] on `grid` [cm⁻¹] for a VECTOR of query points:
column `k` of the result is `σ(grid)` at `pressures[k]` [hPa],
`temperatures[k]` [K] (and, for [`AbscoLUT`](@ref), H₂O broadener `vmr[k]`).
The result is `(length(grid), length(pressures))` on the model's architecture.

The point axis is a flat list: pass one column's layers, or the concatenated
layers of many columns — batching efficiency scales with the total count.

Fused single-kernel implementations exist for [`InterpolationModel`](@ref) and
[`AbscoLUT`](@ref) (`interp = :linear`); anything else falls back to
column-by-column scalar queries with identical results.
"""
function compute_cross_section_profile(model::AbstractCrossSectionModel,
                                       grid::AbstractVector,
                                       pressures::AbstractVector,
                                       temperatures::AbstractVector;
                                       vmr = nothing,
                                       interp::Symbol = :linear)
    n = length(pressures)
    length(temperatures) == n || throw(ArgumentError(
        "compute_cross_section_profile: $(n) pressures vs $(length(temperatures)) temperatures"))
    vmr === nothing || length(vmr) == n || throw(ArgumentError(
        "compute_cross_section_profile: $(n) pressures vs $(length(vmr)) broadener values"))
    # Generic fallback: scalar queries, one column at a time.
    cols = map(1:n) do k
        if vmr === nothing
            compute_cross_section(model, grid, pressures[k], temperatures[k])
        else
            compute_cross_section(model, grid, pressures[k], temperatures[k];
                                  vmr = vmr[k], interp = interp)
        end
    end
    return reduce(hcat, cols)
end

# ---------------------------------------------------------------------------
# shared: host-side bracket precomputation, uploaded once
# ---------------------------------------------------------------------------

# Bracket every query value against ascending `nodes`, clamped exactly like the
# scalar path (`_bracket(nodes, clamp(x, first, last))`). Returns host vectors.
function _bracket_many(nodes, xs::AbstractVector, ::Type{FT}) where {FT}
    n = length(xs)
    lo = Vector{Int32}(undef, n)
    hi = Vector{Int32}(undef, n)
    fr = Vector{FT}(undef, n)
    @inbounds for k in 1:n
        i, i1, f = _bracket(nodes, clamp(FT(xs[k]), first(nodes), last(nodes)))
        lo[k] = Int32(i); hi[k] = Int32(i1); fr[k] = FT(f)
    end
    return lo, hi, fr
end

# In-kernel lower-bound search: the SAME loop as `_lut_resample_kernel!`, so the
# fused kernels inherit its bracket/edge conventions exactly.
@inline function _kernel_lower_bracket(ν, x, n)
    lo = 1
    hi = n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        @inbounds (ν[mid] <= x) ? (lo = mid) : (hi = mid)
    end
    return lo
end

# ---------------------------------------------------------------------------
# InterpolationModel: fused (p, T) bilinear blend + nu resample
# ---------------------------------------------------------------------------

@kernel function _interp_profile_kernel!(out, @Const(grid), @Const(ν), @Const(σ),
                                         @Const(jp), @Const(jp1), @Const(fp),
                                         @Const(kt), @Const(kt1), @Const(ft))
    i, k = @index(Global, NTuple)
    @inbounds begin
        FT = eltype(out)
        x = grid[i]
        n = length(ν)
        if x < ν[1] || x > ν[n]
            out[i, k] = zero(FT)
        else
            # Same 4-term association as the scalar broadcast:
            # ((w00*s00 + w10*s10) + w01*s01) + w11*s11, weights from (fp, ft).
            fpk = fp[k]; ftk = ft[k]
            jpk = jp[k]; jp1k = jp1[k]; ktk = kt[k]; kt1k = kt1[k]
            if n == 1
                out[i, k] = (1 - fpk) * (1 - ftk) * σ[1, jpk, ktk] +
                            fpk * (1 - ftk) * σ[1, jp1k, ktk] +
                            (1 - fpk) * ftk * σ[1, jpk, kt1k] +
                            fpk * ftk * σ[1, jp1k, kt1k]
            else
                lo = _kernel_lower_bracket(ν, x, n)
                f = (x - ν[lo]) / (ν[lo+1] - ν[lo])
                s_lo = (1 - fpk) * (1 - ftk) * σ[lo, jpk, ktk] +
                       fpk * (1 - ftk) * σ[lo, jp1k, ktk] +
                       (1 - fpk) * ftk * σ[lo, jpk, kt1k] +
                       fpk * ftk * σ[lo, jp1k, kt1k]
                s_hi = (1 - fpk) * (1 - ftk) * σ[lo+1, jpk, ktk] +
                       fpk * (1 - ftk) * σ[lo+1, jp1k, ktk] +
                       (1 - fpk) * ftk * σ[lo+1, jpk, kt1k] +
                       fpk * ftk * σ[lo+1, jp1k, kt1k]
                out[i, k] = (1 - f) * s_lo + f * s_hi
            end
        end
    end
end

function compute_cross_section_profile(im::InterpolationModel{FT},
                                       grid::AbstractVector,
                                       pressures::AbstractVector,
                                       temperatures::AbstractVector;
                                       vmr = nothing,
                                       interp::Symbol = :linear) where {FT}
    vmr === nothing || throw(ArgumentError(
        "InterpolationModel tabulates a fixed vmr; do not pass a broadener vector"))
    npt = length(pressures)
    length(temperatures) == npt || throw(ArgumentError(
        "compute_cross_section_profile: $(npt) pressures vs $(length(temperatures)) temperatures"))
    arch = im.architecture
    AT = array_type(arch)
    jp, jp1, fp = _bracket_many(im.p, pressures, FT)
    kt, kt1, ft = _bracket_many(im.T, temperatures, FT)
    Ng = length(grid)
    out = AT(Matrix{FT}(undef, Ng, npt))
    (Ng == 0 || npt == 0) && return out
    gridd = grid === im.ν ? im.ν : AT(collect(FT, grid))
    _interp_profile_kernel!(devi(arch))(out, gridd, im.ν, im.σ,
                                        AT(jp), AT(jp1), AT(fp),
                                        AT(kt), AT(kt1), AT(ft);
                                        ndrange = (Ng, npt))
    synchronize_if_gpu(arch)
    return out
end

# ---------------------------------------------------------------------------
# AbscoLUT: fused (p, T-per-pressure, vmr) blend + nu resample
# ---------------------------------------------------------------------------

@kernel function _absco_profile_kernel!(out, @Const(grid), @Const(ν), @Const(σ),
                                        @Const(jp), @Const(jp1), @Const(fp),
                                        @Const(kv), @Const(kv1), @Const(fv),
                                        @Const(itl), @Const(itl1), @Const(ftl),
                                        @Const(ith), @Const(ith1), @Const(fth))
    i, k = @index(Global, NTuple)
    @inbounds begin
        FT = eltype(out)
        x = grid[i]
        n = length(ν)
        if x < ν[1] || x > ν[n]
            out[i, k] = zero(FT)
        else
            lo = (n == 1) ? 1 : _kernel_lower_bracket(ν, x, n)
            f = (n == 1) ? zero(FT) : (x - ν[lo]) / (ν[lo+1] - ν[lo])
            s_lo = _absco_point(σ, lo, k, jp, jp1, fp, kv, kv1, fv,
                                itl, itl1, ftl, ith, ith1, fth)
            if n == 1
                out[i, k] = s_lo
            else
                s_hi = _absco_point(σ, lo + 1, k, jp, jp1, fp, kv, kv1, fv,
                                    itl, itl1, ftl, ith, ith1, fth)
                out[i, k] = (1 - f) * s_lo + f * s_hi
            end
        end
    end
end

# One table-nu sample for query point k: T-blend inside each bracketing
# pressure's own T axis, vmr-blend, then p-blend — the scalar path's order.
@inline function _absco_point(σ, iν, k, jp, jp1, fp, kv, kv1, fv,
                              itl, itl1, ftl, ith, ith1, fth)
    @inbounds begin
        jpk = jp[k]; jp1k = jp1[k]; fpk = fp[k]
        kvk = kv[k]; kv1k = kv1[k]; fvk = fv[k]
        # lower pressure bracket, its own T bracket
        s0 = (1 - ftl[k]) * σ[iν, kvk, itl[k], jpk] + ftl[k] * σ[iν, kvk, itl1[k], jpk]
        slo = if kvk == kv1k
            s0
        else
            s1 = (1 - ftl[k]) * σ[iν, kv1k, itl[k], jpk] + ftl[k] * σ[iν, kv1k, itl1[k], jpk]
            (1 - fvk) * s0 + fvk * s1
        end
        jpk == jp1k && return slo
        # upper pressure bracket, ITS own T bracket
        t0 = (1 - fth[k]) * σ[iν, kvk, ith[k], jp1k] + fth[k] * σ[iν, kvk, ith1[k], jp1k]
        shi = if kvk == kv1k
            t0
        else
            t1 = (1 - fth[k]) * σ[iν, kv1k, ith[k], jp1k] + fth[k] * σ[iν, kv1k, ith1[k], jp1k]
            (1 - fvk) * t0 + fvk * t1
        end
        return (1 - fpk) * slo + fpk * shi
    end
end

function compute_cross_section_profile(lut::AbscoLUT{FT},
                                       grid::AbstractVector,
                                       pressures::AbstractVector,
                                       temperatures::AbstractVector;
                                       vmr = nothing,
                                       interp::Symbol = :linear) where {FT}
    npt = length(pressures)
    length(temperatures) == npt || throw(ArgumentError(
        "compute_cross_section_profile: $(npt) pressures vs $(length(temperatures)) temperatures"))
    vmr === nothing || length(vmr) == npt || throw(ArgumentError(
        "compute_cross_section_profile: $(npt) pressures vs $(length(vmr)) broadener values"))
    if interp !== :linear
        # :cubic keeps its per-point path (4-node Catmull-Rom per pressure);
        # fall back to the generic column loop rather than silently changing it.
        return invoke(compute_cross_section_profile,
                      Tuple{AbstractCrossSectionModel, AbstractVector,
                            AbstractVector, AbstractVector},
                      lut, grid, pressures, temperatures; vmr = vmr, interp = interp)
    end
    arch = lut.architecture
    AT = array_type(arch)
    jp, jp1, fp = _bracket_many(lut.p, pressures, FT)
    vv = vmr === nothing ? zeros(FT, npt) : vmr
    kv, kv1, fv = _bracket_many(lut.vmr, vv, FT)
    # The ABSCO T axis slides with pressure: bracket T separately inside each
    # bracketing pressure's own axis, exactly as `_absco_at_p` does.
    itl = Vector{Int32}(undef, npt); itl1 = Vector{Int32}(undef, npt); ftl = Vector{FT}(undef, npt)
    ith = Vector{Int32}(undef, npt); ith1 = Vector{Int32}(undef, npt); fth = Vector{FT}(undef, npt)
    @inbounds for k in 1:npt
        Tk = FT(temperatures[k])
        Tlo = @view lut.T[:, jp[k]]
        i, i1, f = _bracket(Tlo, clamp(Tk, first(Tlo), last(Tlo)))
        itl[k] = Int32(i); itl1[k] = Int32(i1); ftl[k] = FT(f)
        Thi = @view lut.T[:, jp1[k]]
        i, i1, f = _bracket(Thi, clamp(Tk, first(Thi), last(Thi)))
        ith[k] = Int32(i); ith1[k] = Int32(i1); fth[k] = FT(f)
    end
    Ng = length(grid)
    out = AT(Matrix{FT}(undef, Ng, npt))
    (Ng == 0 || npt == 0) && return out
    gridd = grid === lut.ν ? lut.ν : AT(collect(FT, grid))
    _absco_profile_kernel!(devi(arch))(out, gridd, lut.ν, lut.σ,
                                       AT(jp), AT(jp1), AT(fp),
                                       AT(kv), AT(kv1), AT(fv),
                                       AT(itl), AT(itl1), AT(ftl),
                                       AT(ith), AT(ith1), AT(fth);
                                       ndrange = (Ng, npt))
    synchronize_if_gpu(arch)
    return out
end
