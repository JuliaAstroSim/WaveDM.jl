# Configs

## Core Configuration Structures

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    SimulationGrid,
    DeviceConfig,
    TimeStepConfig,
    AstroUnitsConfig
]
```

## Kick Step Configuration

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    GravityConfig,
    TidalFieldConfig
]
```

## Initial Conditions Configuration

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    InitialConditionsConfig,
    DensityProfileConfig,
    MassRadiusConfig
]
```

## Visualization Configuration

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    VisualizationConfig,
    VisualizationData
]
```

## Profile Fitting Configuration

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    ProfileFitConfig,
    RCFitConfig,
    BestFitConfig
]
```
