# Introduction

## 1. Wave Dark Matter Background

Wave Dark Matter (WaveDM), also known as Fuzzy Dark Matter or Ultralight Dark Matter, is a theoretical model describing a condensate of extremely light bosonic particles (mass range approximately $10^{-22}$ to $10^{-20}$ eV) at galaxy scales. The de Broglie wavelength of these particles reaches kiloparsec scales, causing quantum interference effects to dominate the dynamics of dark matter halos.

The collective dynamics of WaveDM are described by the Gross-Pitaevskii equation (GPE), which arises from mean-field Hartree-Fock-Bogoliubov theory:

```math
\begin{equation}
    i\hbar \frac{\partial\psi}{\partial t}(\mathbf{r},t) = -\frac{\hbar^2}{2m}\nabla^2 \psi(\mathbf{r},t) + m\Phi_{\mathrm{total}}(\mathbf{r},t)\psi(\mathbf{r},t),
\end{equation}
```

where $\psi$ is the macroscopic wavefunction (order parameter) of the condensate, $m$ is the boson mass, and $\Phi_{\mathrm{total}}$ is the total gravitational potential including contributions from both dark matter and baryonic matter. The gravitational potential is determined by the Poisson equation:

```math
\begin{equation}
    \nabla^2 \Phi_{\mathrm{total}}(\mathbf{r},t) = 4\pi G\left[\rho_{\mathrm{DM}} + \rho_{\mathrm{baryon}}(\mathbf{r},t)\right],
\end{equation}
```

where $\rho_{\mathrm{DM}} = m |\psi|^2$ is the WaveDM density.

WaveDM.jl adopts a wave-based approach, using a split-step Fourier method to solve the Schrödinger-Poisson equation, combining high accuracy with computational efficiency.

## 3. Key Features of WaveDM.jl

### 3.1 High-Performance Numerical Methods

- **Split-step Fourier method**: Accurate solution of wave dynamics without numerical diffusion
- **GPU acceleration**: Support for CUDA acceleration, significantly improving performance
- **Multi-threaded parallelism**: Full utilization of multi-core CPUs
- **Adaptive time-stepping**: Automatic adjustment based on CFL conditions

### 3.2 Flexible Initial Conditions

- Multiple dark matter density profiles: gNFW, Zhao, NFW
- Initial conditions for Milky Way, SPARC galaxies, and dwarf galaxies
- Custom velocity field initialization strategies
- Support for baryonic N-body particles

### 3.3 Realistic Physical Environment

- **Time-dependent tidal fields**: Integration of Milky Way satellite backward trajectories
- **Large Magellanic Cloud perturbations**: Can include gravitational influence of LMC
- **Multiple boundary conditions**: Support for periodic and absorbing boundaries
- **Gravitational potential solvers**: Combination of FFT and N-body methods

### 3.4 Rich Visualization and Analysis

- **Real-time simulation visualization**: Density field, gravitational potential, velocity field
- **Rotation curve analysis**: Automatic generation of rotation curves
- **Density profile fitting**: Support for gNFW and other model fitting
- **Data export**: Support for multiple data formats

## 4. Project Ecosystem

WaveDM.jl integrates with the JuliaAstroSim ecosystem:

- **AstroSimBase**: Basic type definitions
- **PhysicalParticles**: N-body particles and vector algebra
- **PhysicalFFT**: FFT solvers
- **AstroNbodySim**: Gravitational N-body simulations
- **AstroIC**: Initial condition generation
- **AstroPlot**: Visualization tools
- **GalacticDynamics**: Density profiles and orbit calculations
