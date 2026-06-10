#=
The single cross-section kernel. Each workitem owns one grid point, loops over the
active lines whose wing window covers it, and accumulates `S · profile`. `profile`
and `cpf` are singleton type parameters, so the compiler specializes this kernel per
(profile, cpf) with no runtime dispatch — one source, every CPU/CUDA/Metal backend.
=#

@kernel function _crosssection_kernel!(A, @Const(grid), @Const(ν0), @Const(γd), @Const(Γ0),
                                       @Const(Γ2), @Const(Δ0), @Const(Δ2), @Const(νVC),
                                       @Const(η), @Const(Y), @Const(S), @Const(istart), @Const(istop),
                                       N, profile, cpf)
    I  = @index(Global, Linear)
    FT = eltype(A)
    νI = FT(grid[I])
    acc = zero(FT)
    @inbounds for j in 1:N
        if istart[j] ≤ I ≤ istop[j]
            p = (γd = γd[j], Γ0 = Γ0[j], Γ2 = Γ2[j], Δ0 = Δ0[j],
                 Δ2 = Δ2[j], νVC = νVC[j], η = η[j], Y = Y[j])
            acc += S[j] * evaluate(profile, cpf, νI, ν0[j], p)
        end
    end
    @inbounds A[I] += acc
end

"""
    compute_cross_section(model, grid, pressure, temperature;
                          vmr=model.vmr, wavelength_flag=false) -> Vector

Absorption cross-section [cm²/molecule] on `grid` at `pressure` [hPa] and
`temperature` [K]. Result lives on the model's architecture (host `Array` for CPU).

`grid` is wavenumber [cm⁻¹] by default. Set `wavelength_flag=true` to pass `grid` in
**wavelength [nm]** instead: the grid is converted to wavenumber via
`ν[cm⁻¹] = NM_PER_M / λ[nm]` before any line work, so the wing cutoff is always applied
in wavenumber space. The returned `σ` aligns **element-for-element with the input `grid`**
in its original order — a nm-ascending grid (which is wavenumber-descending) comes back in
nm-ascending order, even though the internal computation runs on the sorted wavenumber grid.

`vmr` is the volume mixing ratio of the absorbing gas, used to blend self- and
foreign(air)-broadening: width and shift are `(1-vmr)·foreign + vmr·self`. It defaults to
the model's `vmr`, but can be overridden per call — e.g. to sweep the H₂O cross-section
over humidity without rebuilding the model. `vmr=0` is pure foreign (air) broadening.
"""
function compute_cross_section(model::LineByLineModel{FT}, grid::AbstractVector,
                               pressure::Real, temperature::Real;
                               vmr::Real = model.vmr, wavelength_flag::Bool = false) where {FT}
    arch = model.architecture
    Ng   = length(grid)
    Ng == 0 && return array_type(arch)(zeros(FT, Ng))

    if wavelength_flag
        # Convert nm → cm⁻¹ FIRST (host-side), then run the wavenumber path so the wing
        # cutoff is windowed in cm⁻¹. nm-ascending is cm⁻¹-descending, but `prepare` needs
        # an ascending grid for its binary searches, so we sort, compute, and scatter σ
        # back to the caller's original (nm) order. All ordering work is host-side; only
        # the sorted cm⁻¹ grid is uploaded to the device.
        ν    = FT(NM_PER_M) ./ collect(FT, grid)        # wavenumber [cm⁻¹]
        perm = sortperm(ν)                              # ascending-cm⁻¹ ordering
        σν   = _compute_cross_section(model, ν[perm], pressure, temperature; vmr)
        σh   = Vector{FT}(undef, Ng)                    # inverse-permute on the host
        σh[perm] = Array(σν)                            #   (Array() copies a GPU result back)
        return array_type(arch)(σh)                     # result back on the model's device
    end

    return _compute_cross_section(model, grid, pressure, temperature; vmr)
end

# Wavenumber-grid core. `grid` must be ascending [cm⁻¹] (the LineDatabase / prepare contract).
function _compute_cross_section(model::LineByLineModel{FT}, grid::AbstractVector,
                                pressure::Real, temperature::Real; vmr::Real) where {FT}
    arch = model.architecture
    Ng   = length(grid)
    σ    = array_type(arch)(zeros(FT, Ng))
    Ng == 0 && return σ
    prep = prepare(model, grid, pressure, temperature; vmr)
    if prep.n > 0
        gridd  = array_type(arch)(collect(FT, grid))
        kernel = _crosssection_kernel!(devi(arch))
        kernel(σ, gridd, prep.ν0, prep.γd, prep.Γ0, prep.Γ2, prep.Δ0, prep.Δ2,
               prep.νVC, prep.η, prep.Y, prep.S, prep.istart, prep.istop,
               Int32(prep.n), model.profile, model.cpf; ndrange = Ng)
        synchronize_if_gpu(arch)
    end
    return σ
end
