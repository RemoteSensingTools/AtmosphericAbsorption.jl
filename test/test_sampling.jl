using AtmosphericAbsorption
using Test

@testset "spectral sampling" begin
    @testset "point sampling remains the default ($FT)" for FT in (Float32, Float64)
        model = LineByLineModel(
            oneline_db(FT),
            flatpf(FT);
            profile=Lorentz(),
            wing_cutoff=FT(5),
        )
        grid = collect(FT, 999:FT(0.01):1001)
        default = compute_cross_section(model, grid, FT(1013.25), FT(296))
        explicit = compute_cross_section(
            model,
            grid,
            FT(1013.25),
            FT(296),
            PointSampling(),
        )
        @test default == explicit
        @test eltype(explicit) === FT
        @test isempty(compute_cross_section(
            model,
            FT[],
            FT(1013.25),
            FT(296),
            ConservativeCrossSectionSampling(),
        ))
    end

    @testset "piecewise-linear conservative resampling ($FT)" for FT in (Float32, Float64)
        source_grid = collect(FT, 0:FT(0.01):4)
        linear = @. FT(2) + FT(3) * source_grid
        target_grid = FT[0.5, 1.5, 2.5, 3.5]
        edges = FT[0, 1, 2, 3, 4]
        expected = @. FT(2) + FT(3) * target_grid

        sampled = conservative_resample(
            source_grid,
            linear,
            target_grid;
            cell_edges=edges,
        )
        @test sampled ≈ expected rtol=10eps(FT)
        @test conservative_resample(
            reverse(source_grid),
            reverse(linear),
            reverse(target_grid);
            cell_edges=reverse(edges),
        ) ≈ reverse(expected) rtol=10eps(FT)

        matrix = hcat(linear, FT(2) .* linear)
        matrix_sampled = conservative_resample(
            source_grid,
            matrix,
            target_grid;
            cell_edges=edges,
        )
        @test matrix_sampled[:, 1] ≈ expected rtol=10eps(FT)
        @test matrix_sampled[:, 2] ≈ FT(2) .* expected rtol=10eps(FT)

        @test conservative_resample(
            source_grid,
            linear,
            FT[2];
            cell_edges=FT[1, 3],
        ) ≈ FT[8] rtol=10eps(FT)
        @test size(conservative_resample(source_grid, matrix, FT[])) == (0, 2)
        @test_throws ArgumentError conservative_resample(
            source_grid,
            linear,
            FT[1, 0, 2],
        )
        @test_throws DimensionMismatch conservative_resample(
            source_grid,
            linear,
            target_grid;
            cell_edges=FT[0, 1],
        )
    end

    @testset "cross-section versus transmission conservation ($FT)" for FT in (Float32, Float64)
        source_grid = collect(FT, -0.1:FT(0.0001):0.1)
        narrow_cross_section = @. FT(2) * exp(-(source_grid / FT(0.004))^2)
        target_grid = FT[0]
        edges = FT[-0.1, 0.1]
        column = FT(10)

        mean_cross_section = only(conservative_resample(
            source_grid,
            narrow_cross_section,
            target_grid;
            cell_edges=edges,
        ))
        mean_transmission = only(conservative_resample(
            source_grid,
            exp.(-column .* narrow_cross_section),
            target_grid;
            cell_edges=edges,
        ))
        effective_cross_section = -log(mean_transmission) / column

        # Jensen's inequality: an unresolved opaque spike leaves more light than
        # Beer--Lambert applied to its cell-mean cross-section.
        @test effective_cross_section < mean_cross_section
        @test mean_transmission > exp(-column * mean_cross_section)

        smooth_cross_section = fill(FT(0.02), length(source_grid))
        smooth_mean = only(conservative_resample(
            source_grid,
            smooth_cross_section,
            target_grid;
            cell_edges=edges,
        ))
        smooth_transmission = only(conservative_resample(
            source_grid,
            exp.(-column .* smooth_cross_section),
            target_grid;
            cell_edges=edges,
        ))
        tolerance = FT === Float32 ? FT(2e-4) : FT(100) * eps(FT)
        @test -log(smooth_transmission) / column ≈ smooth_mean rtol=tolerance
    end

    @testset "end-to-end line-by-line conservation ($FT)" for FT in (Float32, Float64)
        model = LineByLineModel(
            oneline_db(FT; S=FT(1e-21), γ_air=FT(0.002)),
            flatpf(FT);
            profile=Lorentz(),
            wing_cutoff=FT(1),
        )
        target_grid = FT[999.9, 1000, 1000.1]
        edges = spectral_cell_edges(target_grid)
        requested_step = minimum(abs, diff(edges)) / FT(100)
        intervals = ceil(Int, (last(edges) - first(edges)) / requested_step)
        fine_grid = collect(range(first(edges), last(edges); length=intervals + 1))
        fine_cross_section = compute_cross_section(
            model,
            fine_grid,
            FT(1013.25),
            FT(296),
        )

        cross_section_sampling = ConservativeCrossSectionSampling(refinement=100)
        conserved_cross_section = compute_cross_section(
            model,
            target_grid,
            FT(1013.25),
            FT(296),
            cross_section_sampling,
        )
        reference_cross_section = conservative_resample(
            fine_grid,
            fine_cross_section,
            target_grid;
            cell_edges=edges,
        )
        @test conserved_cross_section ≈ reference_cross_section rtol=50eps(FT)
        @test eltype(conserved_cross_section) === FT

        column = FT(1e21)
        transmission_sampling = ConservativeTransmissionSampling(
            column;
            refinement=100,
        )
        effective_cross_section = compute_cross_section(
            model,
            target_grid,
            FT(1013.25),
            FT(296),
            transmission_sampling,
        )
        reference_transmission = conservative_resample(
            fine_grid,
            exp.(-column .* fine_cross_section),
            target_grid;
            cell_edges=edges,
        )
        @test exp.(-column .* effective_cross_section) ≈
              reference_transmission rtol=100eps(FT)
        @test all(effective_cross_section .<= conserved_cross_section)
        @test eltype(effective_cross_section) === FT

        nm_per_m = FT(AtmosphericAbsorption.Constants.NM_PER_M)
        wavelength_model = LineByLineModel(
            oneline_db(FT; ν0=nm_per_m / FT(770.4)),
            flatpf(FT);
            profile=Lorentz(),
            wing_cutoff=FT(5),
        )
        wavelength_grid = FT[770.2, 770.4, 770.6]
        wavelength_cross_section = compute_cross_section(
            wavelength_model,
            wavelength_grid,
            FT(1013.25),
            FT(296),
            ConservativeCrossSectionSampling(fine_step=FT(0.001));
            wavelength_flag=true,
        )
        @test argmax(wavelength_cross_section) == 2
        @test eltype(wavelength_cross_section) === FT
    end

    @testset "sampling configuration validation" begin
        @test_throws ArgumentError ConservativeCrossSectionSampling(refinement=0)
        @test_throws ArgumentError ConservativeCrossSectionSampling(fine_step=0)
        @test_throws ArgumentError ConservativeTransmissionSampling(0.0)
        @test_throws ArgumentError ConservativeTransmissionSampling(Inf)
    end
end
