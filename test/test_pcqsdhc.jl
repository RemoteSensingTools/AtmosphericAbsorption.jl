using AtmosphericAbsorption.LineShapes: pcqsdhc, voigt, HumlicekWeideman32, ErfcxCPF
using DelimitedFiles: readdlm
using JET
# GOLDEN (the reference_data artifact's golden dir) comes from runtests.jl.

# name => (Γ0, Γ2, Δ0, Δ2, νVC, η, Ylm); goldens from generate_pcqsdhc_golden.py (HAPI).
# The "htlm" case carries first-order line mixing; golden col2 = real(LS) + Ylm·imag(LS).
const PCQSDHC_CASES = Dict(
    "voigt" => (0.03, 0.0, -0.01, 0.0, 0.0, 0.0, 0.0),
    "sdv"   => (0.03, 0.006, -0.01, 0.002, 0.0, 0.0, 0.0),
    "ht"    => (0.03, 0.006, -0.01, 0.002, 0.012, 0.3, 0.0),
    "htlm"  => (0.03, 0.006, -0.01, 0.002, 0.012, 0.3, 0.4),
)

@testset "pcqsdhc (HT/SDV/line-mixing) vs HAPI golden" begin
    ν0, γd = 1000.0, 0.02
    @testset "$name ($FT)" for name in keys(PCQSDHC_CASES), FT in (Float32, Float64)
        Γ0, Γ2, Δ0, Δ2, νVC, η, Ylm = PCQSDHC_CASES[name]
        d = readdlm(joinpath(GOLDEN, "pcqsdhc_$(name).txt"), comments = true, comment_char = '#')
        ν, re_ref = d[:, 1], d[:, 2]
        re = [let LS = pcqsdhc(HumlicekWeideman32(), FT(ν0), FT(γd), FT(Γ0), FT(Γ2),
                               FT(Δ0), FT(Δ2), FT(νVC), FT(η), FT(νi))
                  real(LS) + FT(Ylm) * imag(LS)            # first-order line-mixing fold
              end for νi in ν]
        rel = maximum(abs.(re .- re_ref)) / maximum(abs.(re_ref))
        @test rel < (FT === Float64 ? 1e-4 : 1e-3)
    end

    @testset "HT reduces to Voigt when Γ2=Δ2=νVC=η=0 ($FT)" for FT in (Float32, Float64)
        γd, Γ0 = FT(0.02), FT(0.03)
        y = sqrt(log(FT(2))) * Γ0 / γd
        cpf = ErfcxCPF()
        for νi in FT.(999.6:0.01:1000.4)
            ht = real(pcqsdhc(cpf, FT(1000), γd, Γ0, FT(0), FT(0), FT(0), FT(0), FT(0), νi))
            @test isapprox(ht, voigt(cpf, νi - FT(1000), γd, y); rtol = FT(1e-5))
        end
    end

    @testset "type stability ($FT)" for FT in (Float32, Float64)
        z = @inferred pcqsdhc(HumlicekWeideman32(), FT(1000), FT(0.02), FT(0.03),
                              FT(0.006), FT(-0.01), FT(0.002), FT(0.012), FT(0.3), FT(1000.1))
        @test z isa Complex{FT}
    end
end
