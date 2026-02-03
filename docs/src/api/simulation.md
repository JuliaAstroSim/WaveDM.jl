# Simulation

```@meta
CurrentModule = WaveDM
```

## Main Simulation Functions

### SPE3D_waveDM

```@docs
WaveDM.SPE3D_waveDM
```

### simulate_waveDM

```@docs
WaveDM.simulate_waveDM
```

## Simulation Overview

WaveDM.jl provides a high-level interface for simulating Wave Dark Matter dynamics using the Schrödinger-Poisson equation. The main entry points are the `SPE3D_waveDM` and `simulate_waveDM` functions, which handle the complete simulation workflow from initial condition setup to result visualization.

### Simulation Workflow

1. **Grid Setup**: Create a uniform Cartesian grid with specified dimensions and resolution
2. **Initial Conditions**: Generate initial wavefunction based on density profile
3. **Velocity Field**: Create initial velocity field consistent with the gravitational potential
4. **Configuration**: Set up simulation parameters, boundary conditions, and visualization options
5. **Evolution**: Run the split-step Fourier method to evolve the wavefunction
6. **Visualization**: Generate real-time or post-simulation plots of density, potential, and other quantities
7. **Analysis**: Compute simulation statistics and fit results to target profiles

### Key Simulation Parameters

#### Grid and Resolution

- `Xmax`, `Ymax`, `Zmax`: Simulation box dimensions (in dimensionless units)
- `Nx`, `Ny`, `Nz`: Number of grid points in each dimension
- `Tmax`: Total simulation time (in dimensionless units)
- `Nt`: Number of time steps

#### Initial Conditions

- `IC`: Initial wavefunction or distribution function
- `V`: External potential function
- `κ`: Self-interaction parameter
- `baryon`: Baryonic density contribution

#### Physical Models

- `baryon_mode`: How baryons are treated (:mesh, :particles_static, :ignored)
- `boundary`: Boundary conditions (Periodic() or Vacuum())
- `GravitySolver`: Gravitational solver type (Tree() or DirectSummation())

#### Performance Options

- `gpu`: Enable GPU acceleration
- `distributed_memory`: Enable distributed memory parallelism
- `autoset_timestep`: Automatically set timestep based on stability criteria

#### Visualization

- `title`: Simulation title for output files
- `outputdir`: Directory for saving output files
- `Realtime`: Enable real-time visualization
- `StepsBetweenSnapshots`: Frequency of visualization updates

### Example Usage

#### Basic Simulation

```julia
using WaveDM

# Run a simple simulation of an isolated halo
simulate_waveDM(
    model = :dwarf,
    title = "QuickTest",
    Xmax = 2.0,
    Tmax = 0.001,
    Nt = 3,
    Nx = Ny = Nz = 64,
    Np = 100,
    gpu = false
)
```

#### Advanced Simulation with Custom Parameters

```julia
using WaveDM

# Run a more detailed simulation with custom parameters
SPE3D_waveDM(
    Xmax = 5.0,
    Ymax = 5.0,
    Zmax = 5.0,
    Tmax = 1.0,
    Nx = 128,
    Ny = 128,
    Nz = 128,
    Nt = 1000,
    autoset_timestep = true,
    IC = (x,y,z)->exp(-x^2-y^2-z^2),  # Gaussian initial condition
    V = (x,y,z,ψ)->0.0,  # No external potential
    baryon_mode = :mesh,  # Baryons on mesh
    gpu = true,  # Enable GPU acceleration
    Realtime = true,  # Show real-time plots
    StepsBetweenSnapshots = 10,  # Update plots every 10 steps
    title = "DetailedWaveDMSim",
    outputdir = joinpath(pwd(), "output")
)
```

### Output Files

The simulation generates several types of output files:

1. **Visualization Files**:
   - PNG images of density slices and profiles
   - MP4 video of simulation evolution

2. **Data Files**:
   - JLD2 files containing wavefunction and potential data
   - CSV files with simulation statistics
   - Configuration files for reproducibility

3. **Analysis Files**:
   - Fitted profiles and rotation curves
   - Best-fit results for dwarf galaxy simulations

### Performance Considerations

- **GPU Acceleration**: Enabling `gpu=true` can provide significant speedups for large simulations
- **Grid Size**: Powers of two (64, 128, 256, etc.) are optimal for FFT performance
- **Time Step**: Smaller time steps improve accuracy but increase computational cost
- **Visualization**: Disabling `Realtime` can speed up simulations by avoiding plotting overhead

## Related Modules

- [Configs](configs.md): Simulation configuration structures
- [KDK](KDK.md): Split-step Fourier method implementation