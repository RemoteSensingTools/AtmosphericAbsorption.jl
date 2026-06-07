# API reference

```@meta
CurrentModule = AtmosphericAbsorption
```

## Line lists & data sources

```@docs
LineDatabase
SourceMetadata
load_lines
source_metadata
HitranPort
fetch_hitran
load_hitran
ExoMolPort
activate_hitran!
fetch_hitran_nonvoigt
load_hitran_nonvoigt
```

## Species notation

Name molecules and isotopologues with the backend-agnostic notation (see the
[Species & isotopologues](isotopologues.md) table); the mapping is overridable.

```@docs
molecules
molecule_number
molecule_symbol
isotopologue
register_molecule!
register_isotopologue!
```

## Cross-section model

```@docs
LineByLineModel
compute_cross_section
```

## Interpolation tables (precomputed LUTs)

```@docs
InterpolationModel
build_interpolation_model
save_interpolation_model
```

## Tabulated cross-sections (HITRAN `.xsc`)

```@docs
TabulatedCrossSection
XscBand
read_xsc
load_xsc
fetch_hitran_xsc
```

## Line shapes

```@docs
AtmosphericAbsorption.LineShapes.Doppler
AtmosphericAbsorption.LineShapes.Lorentz
AtmosphericAbsorption.LineShapes.Voigt
AtmosphericAbsorption.LineShapes.SpeedDependentVoigt
AtmosphericAbsorption.LineShapes.Rautian
AtmosphericAbsorption.LineShapes.SpeedDependentRautian
AtmosphericAbsorption.LineShapes.HartmannTran
AtmosphericAbsorption.LineShapes.doppler
AtmosphericAbsorption.LineShapes.lorentz
AtmosphericAbsorption.LineShapes.voigt
AtmosphericAbsorption.LineShapes.HumlicekWeideman32
AtmosphericAbsorption.LineShapes.ErfcxCPF
```

## Partition functions (advanced)

Normally the partition function rides on the `LineDatabase` (set by the port) and the
model picks it up automatically — you only reach for these to pin a specific edition.

```@docs
partition_function
AbstractPartitionFunction
TIPS2021PF
TIPS2017PF
TabulatedPF
Q_ratio
pf_name
```

## Continuum

```@docs
load_mtckd
build_mtckd_band
h2o_continuum
h2o_continuum!
load_cia
parse_cia_file
cia_cross_section
cia_cross_section!
```

## Compute backends

```@docs
AtmosphericAbsorption.Architectures.CPU
AtmosphericAbsorption.Architectures.GPU
AtmosphericAbsorption.Architectures.MetalGPU
AtmosphericAbsorption.Architectures.default_architecture
```
