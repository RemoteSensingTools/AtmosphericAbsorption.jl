"""Supertype for cross-section models (line-by-line, interpolated, …)."""
abstract type AbstractCrossSectionModel end

"""
    LineByLineModel(lines, partition; profile=Voigt(), cpf=HumlicekWeideman32(),
                    wing_cutoff=40.0, vmr=0.0, architecture=default_architecture())

A line-by-line cross-section model: a `LineDatabase`, its partition function, the
line `profile` and complex-probability strategy `cpf`, plus the wing cutoff [cm⁻¹],
broadening `vmr`, and compute `architecture`. Singleton `profile`/`cpf` types let the
kernel specialize at compile time.
"""
struct LineByLineModel{FT,L<:LineDatabase{FT},P<:AbstractLineProfile,
                       C<:AbstractCPF,PF<:AbstractPartitionFunction,
                       A<:AbstractArchitecture} <: AbstractCrossSectionModel
    lines::L
    partition::PF
    profile::P
    cpf::C
    wing_cutoff::FT
    vmr::FT
    architecture::A
end

function LineByLineModel(lines::LineDatabase{FT}, partition::AbstractPartitionFunction;
                         profile::AbstractLineProfile = Voigt(),
                         cpf::AbstractCPF = HumlicekWeideman32(),
                         wing_cutoff = 40, vmr = 0,
                         architecture::AbstractArchitecture = default_architecture()) where {FT}
    return LineByLineModel(lines, partition, profile, cpf,
                           FT(wing_cutoff), FT(vmr), architecture)
end
