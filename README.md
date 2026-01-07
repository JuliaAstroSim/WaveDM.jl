# WaveDM.jl

[![codecov](https://codecov.io/gh/JuliaAstroSim/WaveDM.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/WaveDM.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastrosim.github.io/WaveDM.jl/)

WaveDM.jl is an open-source Julia package for simulating ultralight wave dark matter (WaveDM) dynamics at galaxy scales. The code solves the time-dependent Schrödinger-Poisson equation (SPE) using a split-step Fourier method, supporting various initial conditions, boundary conditions, and external gravitational perturbations.

## Project Background

Wave dark matter (also known as fuzzy dark matter or ultralight dark matter) is a theoretical model proposing that dark matter consists of extremely light bosonic particles with masses in the range $10^{-22}$ to $30$ eV. At these masses ($\sim 10^{-22}$ eV), the de Broglie wavelength of dark matter particles reaches kiloparsec scales, causing quantum interference effects to dominate the dynamics of dark matter halos.

WaveDM.jl aims to provide a high-performance, user-friendly framework for simulating WaveDM dynamics, addressing several limitations of existing codes:
- Flexible parallelization strategies (GPU, multi-threading, distributed computing)
- Rich initial condition generators for diverse galactic systems
- Support for time-dependent tidal fields and external perturbations
- Multiple boundary conditions (periodic, absorbing)
- Real-time visualization and interactive analysis tools

## Installation

1. Install Julia from the [official website](https://julialang.org/downloads/)
2. Open Julia and enter the package manager by pressing `]`
3. Install WaveDM.jl and other dependencies:
   ```julia
   pkg> add WaveDM
   ```

### For Development

```bash
git clone https://github.com/JuliaAstroSim/WaveDM.jl.git
cd WaveDM.jl
julia --project -e "using Pkg; Pkg.develop(Pkg.PackageSpec(path='.'))"
```

## Quick Start

### Schrödinger-Poisson Equation

WaveDM.jl solves the time-dependent Schrödinger-Poisson equation:

```math
\begin{align}
    i\hbar \frac{\partial\psi}{\partial t} &= -\frac{\hbar^2}{2m} \nabla^2 \psi + m\Phi \psi, \\
    \nabla^2 \Phi &= 4\pi G m |\psi|^2,
\end{align}
```

where $\psi$ is the macroscopic wavefunction, $m$ is the boson mass, and $\Phi$ is the gravitational potential.

### Split-step Fourier Method

The code uses a second-order split-step Fourier method with a "kick-drift-kick" leapfrog scheme:

```math
\psi(\mathbf{r}, t + \Delta t) \approx e^{i\frac{\Delta t}{2}\hat{H}_p} e^{i\Delta t \hat{H}_k} e^{i\frac{\Delta t}{2}\hat{H}_p} \psi(\mathbf{r}, t)
```

where $\hat{H}_k = -\frac{1}{2}\nabla^2$ is the kinetic energy operator and $\hat{H}_p = \Phi$ is the potential energy operator.

### Gravitational Potential Solvers

- **FFT-based solver** for WaveDM density field
- **Tree-based solver** for baryonic N-body particles
- **Direct summation** for small particle counts

## Usage Examples

### Isolated Halo Simulation

### Baryonic Component Integration

### Documentation

Comprehensive documentation is available at [https://juliaastrosim.github.io/WaveDM.jl/](https://juliaastrosim.github.io/WaveDM.jl/).

## Performance Considerations

- **GPU Acceleration**: Enable `gpu=true` 
- **Grid Size**: Powers of two (64, 128, 256, etc.) are optimal for FFT performance
- **Visualization**: Adjust `StepsBetweenSnapshots` to balance visualization and performance
- **Memory Usage**: Monitor GPU memory usage for large simulations, reduce grid size if needed



## Citation

If you use WaveDM.jl in your research, please cite:

```
@article{WaveDMjl,
  title={WaveDM.jl: A High-Performance Schrödinger-Poisson Solver for Wave Dark Matter Dynamics at Galaxy Scale},
  author={Meng, Run-Yu and Dong, Xiao-Bo},
  journal={Research in Astronomy and Astrophysics},
  year={2026},
  volume={X},
  pages={000--000}
}
```

## License



## Acknowledgments

WaveDM.jl development was supported by:
- Yunnan Astronomical Observatory, National Astronomical Observatories, Chinese Academy of Sciences
- National Natural Science Foundation of China (NSFC)

## Related Packages

WaveDM.jl is part of the JuliaAstroSim ecosystem:
- [AstroSimBase.jl](https://github.com/JuliaAstroSim/AstroSimBase.jl) - Basic types and interfaces
- [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl) - N-body particles and vector algebra
- [AstroNbodySim.jl](https://github.com/JuliaAstroSim/AstroNbodySim.jl) - Gravitational N-body simulations
- [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl) - Initial condition generation
- [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) - Visualization tools

## Contact

For questions, issues, or feature requests, please:
- Open an issue on [GitHub](https://github.com/JuliaAstroSim/WaveDM.jl/issues)

## Roadmap

