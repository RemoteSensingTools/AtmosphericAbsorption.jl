#=
The single cross-section kernel. Each workitem owns one grid point, loops over the
active lines whose wing window covers it, and accumulates `S · profile`. `profile`
and `cpf` are singleton type parameters, so the compiler specializes this kernel per
(profile, cpf) with no runtime dispatch — one source, every CPU/CUDA/Metal backend.
=#

@kernel function _crosssection_kernel!(A, @Const(grid), @Const(ν0), @Const(γd), @Const(Γ0),
                                       @Const(Γ2), @Const(Δ0), @Const(Δ2), @Const(νVC),
                                       @Const(η), @Const(S), @Const(istart), @Const(istop),
                                       N, profile, cpf)
    I  = @index(Global, Linear)
    FT = eltype(A)
    νI = FT(grid[I])
    acc = zero(FT)
    @inbounds for j in 1:N
        if istart[j] ≤ I ≤ istop[j]
            p = (γd = γd[j], Γ0 = Γ0[j], Γ2 = Γ2[j], Δ0 = Δ0[j],
                 Δ2 = Δ2[j], νVC = νVC[j], η = η[j])
            acc += S[j] * evaluate(profile, cpf, νI, ν0[j], p)
        end
    end
    @inbounds A[I] += acc
end

"""
    compute_cross_section(model, grid, pressure, temperature) -> Vector

Absorption cross-section [cm²/molecule] on `grid` [cm⁻¹] at `pressure` [hPa] and
`temperature` [K]. Result lives on the model's architecture (host `Array` for CPU).
"""
function compute_cross_section(model::LineByLineModel{FT}, grid::AbstractVector,
                               pressure::Real, temperature::Real) where {FT}
    arch = model.architecture
    Ng   = length(grid)
    σ    = array_type(arch)(zeros(FT, Ng))
    Ng == 0 && return σ
    prep = prepare(model, grid, pressure, temperature)
    if prep.n > 0
        gridd  = array_type(arch)(collect(FT, grid))
        kernel = _crosssection_kernel!(devi(arch))
        kernel(σ, gridd, prep.ν0, prep.γd, prep.Γ0, prep.Γ2, prep.Δ0, prep.Δ2,
               prep.νVC, prep.η, prep.S, prep.istart, prep.istop,
               Int32(prep.n), model.profile, model.cpf; ndrange = Ng)
        synchronize_if_gpu(arch)
    end
    return σ
end
