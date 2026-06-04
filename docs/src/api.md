# API reference

```@meta
CurrentModule = AtmosphericAbsorption
```

## Line lists & data sources

```@docs
LineDatabase
SourceMetadata
load_lines
partition_function
source_metadata
HitranPort
fetch_hitran
ExoMolPort
activate_hitran!
fetch_hitran_nonvoigt
load_hitran_nonvoigt
```

## Cross-section model

```@docs
LineByLineModel
compute_cross_section
```

## Line shapes

```@docs
AtmosphericAbsorption.LineShapes.Doppler
AtmosphericAbsorption.LineShapes.Lorentz
AtmosphericAbsorption.LineShapes.Voigt
AtmosphericAbsorption.LineShapes.SpeedDependentVoigt
AtmosphericAbsorption.LineShapes.HartmannTran
AtmosphericAbsorption.LineShapes.doppler
AtmosphericAbsorption.LineShapes.lorentz
AtmosphericAbsorption.LineShapes.voigt
AtmosphericAbsorption.LineShapes.HumlicekWeideman32
AtmosphericAbsorption.LineShapes.ErfcxCPF
```

## Partition functions

```@docs
AbstractPartitionFunction
TIPS2017PF
TabulatedPF
Q_ratio
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
