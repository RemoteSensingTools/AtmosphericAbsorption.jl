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

function LineByLineModel(lines::LineDatabase{FT}, partition::AbstractPartitionFunction = lines.partition;
                         profile::AbstractLineProfile = Voigt(),
                         cpf::AbstractCPF = HumlicekWeideman32(),
                         wing_cutoff = 40, vmr = 0,
                         architecture::AbstractArchitecture = default_architecture()) where {FT}
    return LineByLineModel(lines, partition, profile, cpf,
                           FT(wing_cutoff), FT(vmr), architecture)
end

function Base.show(io::IO, m::LineByLineModel{FT}) where {FT}
    print(io, "LineByLineModel{", FT, "}(", length(m.lines), " lines, ",
          nameof(typeof(m.profile)), ", ", pf_name(m.partition), ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::LineByLineModel{FT}) where {FT}
    print(io, "LineByLineModel{", FT, "}")
    print(io, "\n  lines:        ", length(m.lines), " transitions — ", join(molecules(m.lines), ", "))
    print(io, "\n  source:       ", m.lines.meta.source)
    print(io, "\n  partition:    ", pf_name(m.partition))
    print(io, "\n  profile:      ", nameof(typeof(m.profile)))
    print(io, "\n  cpf:          ", nameof(typeof(m.cpf)))
    print(io, "\n  wing cutoff:  ", m.wing_cutoff, " cm⁻¹")
    print(io, "\n  vmr:          ", m.vmr)
    print(io, "\n  architecture: ", nameof(typeof(m.architecture)))
end
