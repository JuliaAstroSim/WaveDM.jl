# WaveDM.jl

WaveDM.jl is an open-source Julia package for simulating ultralight wave dark matter (WaveDM) dynamics at galaxy scales. The code solves the time-dependent Schrödinger-Poisson equation (SPE) using a split-step Fourier method, supporting various initial conditions and boundary conditions, as well as tidal interactions and external gravitational perturbations.

## Core Features

- Solves the 3D time-dependent Schrödinger-Poisson equation
- Supports multiple boundary conditions (periodic, absorbing)
- Flexible initial condition generators
- Integrates N-body gravitational solvers
- Supports time-dependent tidal fields
- Real-time visualization and post-processing tools
- GPU acceleration support

## Quick Start

```julia
using WaveDM

# Run minimal example
test_MW_MOND(;
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

## Key Features

1. **High-Performance Numerical Methods**
   - Split-step Fourier method for accurate wave dynamics
   - GPU acceleration support
   - Multi-threaded parallel computing

2. **Flexible Initial Conditions**
   - Supports multiple dark matter density profiles (gNFW, Zhao, NFW)
   - Generates initial conditions for Milky Way, SPARC galaxies, and dwarfs
   - Custom velocity field initialization

3. **Realistic Physical Environment**
   - Time-dependent tidal fields
   - Backward trajectory for Milky Way satellites
   - Can include gravitational perturbations from the Large Magellanic Cloud

4. **Rich Visualization**
   - Real-time simulation visualization
   - Density field, gravitational potential, velocity field plots
   - Rotation curve and density profile analysis

## Documentation Navigation

- [Introduction](introduction.md) - Wave dark matter background and project motivation
- [Installation](installation.md) - Installing WaveDM.jl and dependencies
- [Algorithms](algorithms.md) - Schrödinger-Poisson equation and numerical methods
- [API Reference](api/configs.md) - Detailed API documentation
- [Examples](examples.md) - Complete simulation examples
- [References](reference.md) - Related literature

## Citation

If you use WaveDM.jl in your research, please cite the following paper:

```
@article{WaveDM.jl,
  title={WaveDM.jl: A High-Performance Schr{"o}dinger-Poisson Solver for Wave Dark Matter Dynamics at Galaxy Scale},
  author={Meng, Run-Yu and Dong, Xiao-Bo},
  journal={Research in Astronomy and Astrophysics},
  year={2026},
  volume={X},
  pages={000--000}
}
```

## License

WaveDM.jl is released under the MIT License. See the [LICENSE](https://github.com/JuliaAstroSim/WaveDM.jl/blob/main/LICENSE) file for details.
