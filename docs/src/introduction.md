# Introduction

## 1. Wave Dark Matter Background

Wave (or "fuzzy") cold dark matter describes dark matter as an ultralight bosonic field with particle masses in the range $m_a \sim 10^{-23}$–$10^{-20}$ eV [Hu et al. 2000; Hui 2021].
For such small masses, the de Broglie wavelength

```math
\lambda_{\rm dB} = \frac{h}{m_a\,v}
```

reaches astrophysical scales inside galactic halos. For a typical halo velocity $v \sim 100$ km s$^{-1}$ and $m_a = 10^{-22}$ eV, $\lambda_{\rm dB} \approx 1.2$ kpc — comparable to the size of galactic substructure. Wave CDM therefore exhibits phenomena absent in particle CDM: solitonic cores supported by quantum pressure, granule structures from wave interference, halo fluctuations, and the suppression of small-scale power below $\lambda_{\rm dB}$.

### 1.1 Bose–Einstein condensation and the classical-field limit

For ultralight bosons the critical temperature for Bose–Einstein condensation is extremely high,
$T_c \sim 10^{35}$ K,
while the cosmic microwave background temperature is $\sim 3$ K and the effective kinetic temperature of the bosons is only $T_{\rm eff} \sim 10^{-25}$ K for $m_a \sim 10^{-22}$ eV. The system is therefore in the effective zero-temperature regime.

In this regime the dynamics of a highly Bose-degenerate bosonic system is well described by a single complex-valued macroscopic wavefunction $\psi(\mathbf r, t)$, the *order parameter* of the condensate. Once equilibrium is reached, almost all bosons reside in the condensate; the ground-state component of the condensate is the *soliton* familiar from the wave-CDM literature, while the excited component carries collective modes (phonons, vortices, localized modes).

Crucially, the classical-field description is not limited to systems that have already established global phase coherence. Even far from equilibrium, when the bosons' de Broglie wavelengths overlap strongly and many low-energy modes are highly occupied, the system can still be evolved with the classical-field GPE — interpreting $\psi(\mathbf r,t)$ as a classical matter-wave field. This is the regime relevant for wave-CDM halos.

### 1.2 The Schrödinger–Poisson system

The evolution of $\psi$ is governed by the Gross–Pitaevskii equation (GPE), coupled self-consistently with the Poisson equation of gravity in the case of DM halos:

```math
\begin{aligned}
i\hbar\,\partial_t \psi &= -\frac{\hbar^2}{2m}\nabla^2\psi + m\,\Phi_{\rm total}\,\psi, \\
\nabla^2 \Phi_{\rm total} &= 4\pi G\bigl[\rho_{\rm DM} + \rho_{\rm baryon}\bigr],
\end{aligned}
```

with $\rho_{\rm DM} = m|\psi|^2$ and $\Phi_{\rm total}$ the total gravitational potential.

Applying the Madelung transformation $\psi = \sqrt{\rho/m}\,e^{iS/\hbar}$ yields a quantum hydrodynamic form: a continuity equation and a quantum Bernoulli equation, where the *quantum potential*

```math
Q = -\frac{\hbar^2}{2m}\,\frac{\nabla^2\sqrt{\rho}}{\sqrt{\rho}}
```

encodes the Heisenberg uncertainty principle. The corresponding *quantum pressure* $-(1/m)\nabla Q$ enters the quantum Euler equation.

A density-weighted volume integral of $Q$ defines the **quantum (gradient) energy**

```math
\mathcal{Q} \;\equiv\; \int \rho Q\,d^3r \;=\; \frac{\hbar^2}{2m}\int |\nabla\sqrt{\rho}|^2 d^3r,
```

where the surface term vanishes for any realistic NFW-like halo. The total energy $\mathcal{E} = \mathcal{K} + \mathcal{V} + \mathcal{Q}$ is conserved for an isolated, fully self-gravitating halo.

### 1.3 Dimensionless units

Following Edwards et al. (2018) and Glennon et al. (2021), WaveDM.jl adopts the characteristic units

```math
\begin{aligned}
\mathcal{L} &= \left(\frac{8\pi\hbar^2}{3m^2 H_0^2\Omega_{m_0}}\right)^{1/4}
            \approx 121\left(\frac{10^{-23}\,\mathrm{eV}}{m}\right)^{1/2}\,\mathrm{kpc}, \\
\mathcal{T} &= \left(\frac{8\pi}{3H_0^2\Omega_{m_0}}\right)^{1/2}
            \approx 75.5\,\mathrm{Gyr}, \\
\mathcal{M} &= \frac{1}{G}\left(\frac{8\pi}{3H_0^2\Omega_{m_0}}\right)^{-1/4}
            \left(\frac{\hbar}{m}\right)^{3/2}
            \approx 7\times10^7\left(\frac{10^{-23}\,\mathrm{eV}}{m}\right)^{3/2}\,M_\odot.
\end{aligned}
```

In these units the SPE takes the compact form

```math
\begin{aligned}
i\,\partial_t \psi &= -\tfrac{1}{2}\nabla^2 \psi + \Phi_{\rm total}\,\psi, \\
\nabla^2 \Phi_{\rm total} &= 4\pi\bigl(\rho_{\rm DM} + \rho_{\rm baryon}\bigr),
\end{aligned}
```

with $\psi = \sqrt{\rho}\,e^{i\theta}$ and $\mathbf v = \nabla\theta$.
The self-interaction term $\kappa|\psi|^2\psi$, with $\kappa = 4\pi\hbar\,a_s/(\mathcal{T} m^2)$, can be added to model repulsive ($\kappa > 0$) or attractive ($\kappa < 0$) boson-boson contact interactions. WaveDM.jl supports $\kappa = 0$ (non-interacting), $\kappa > 0$, and $\kappa < 0$ regimes.

### 1.4 A generalized nonlinear Schrödinger interface

For cross-disciplinary use, WaveDM.jl solves the **generalized nonlinear Schrödinger equation (GNLSE)**

```math
i\,\partial_t \psi(\mathbf r, t) = -\tfrac{1}{2}\nabla^2\psi + V(\mathbf r, t, \psi)\,\psi,
```

where $V$ is an arbitrary user-supplied potential. Two classes of $V$ are supported natively:

- **Local interactions** — $V = F(|\psi|)$ depending pointwise on the field amplitude. Example: the Kerr nonlinearity $V = \kappa|\psi|^2$ for self-interacting wave CDM and instantaneous nonlinear optical response.
- **Nonlocal interactions** — $V = \int U(\mathbf r - \mathbf r')\,|\psi(\mathbf r')|^2 d^3 r'$, or equivalently the solution of an auxiliary field equation such as the Poisson equation. Example: the self-gravitational potential in wave CDM, dipole–dipole interactions in atomic condensates, thermal-optical nonlinearities in nonlinear optics.

This interface is the foundation of WaveDM.jl's cross-disciplinary applicability: astrophysics, nonlinear optics, and cold-atom physics all reduce to choosing an appropriate $V$ — the same code base, the same solver, the same parallel infrastructure.

## 2. Design Philosophy

The primary design goal of WaveDM.jl is to **make galaxy-scale simulations accessible to a broad range of astrophysicists** who are not necessarily computational specialists, while remaining a community-oriented open-source project that **advanced users and developers can extend** without mastering the entire codebase.

### 2.1 Adaptability for general users

Wave dark matter simulations typically demand expertise across numerical methods for the Schrödinger equation, parallel computing, initial conditions, N-body gravity, and astrophysical analysis. WaveDM.jl unifies these components behind a small set of high-level keyword arguments.

A typical user might want to:

- vary the boson mass $m_a$ and self-interaction strength $\kappa$ to explore different wave CDM scenarios;
- enable or disable tidal perturbations from the LMC or a host halo;
- change numerical resolution, the initial halo, or the velocity-field configuration.

In every case the user adjusts a few keywords — they do **not** edit the underlying implementation. This is enabled by:

- the **generalized SPE** formulation $V(\mathbf r, t, \psi)\psi$, which exposes a single pluggable potential interface;
- the **galaxy-simulation toolbox** (Section 3) with ready-to-use initial conditions, trajectory lookback, tidal fields, and real-time visualization.

### 2.2 Extensibility for advanced users

WaveDM.jl is built on the [JuliaAstroSim](https://github.com/JuliaAstroSim) ecosystem of reusable components. New physical models — additional baryonic feedback mechanisms, SPH for gas dynamics, etc. — can be added by implementing a small set of well-defined interfaces. The shared infrastructure includes:

- [`PhysicalParticles.jl`](https://github.com/JuliaAstroSim/PhysicalParticles.jl) — fundamental data structures and vector algebra;
- [`AstroIO.jl`](https://github.com/JuliaAstroSim/AstroIO.jl) — file I/O;
- [`AstroIC.jl`](https://github.com/JuliaAstroSim/AstroIC.jl) — initial conditions;
- [`AstroPlot.jl`](https://github.com/JuliaAstroSim/AstroPlot.jl) — analysis and visualization;
- [`PhysicalTrees.jl`](https://github.com/JuliaAstroSim/PhysicalTrees.jl) and [`PhysicalMeshes.jl`](https://github.com/JuliaAstroSim/PhysicalMeshes.jl) — tree and mesh gravity solvers;
- [`PhysicalFDM.jl`](https://github.com/JuliaAstroSim/PhysicalFDM.jl) and [`PhysicalFFT.jl`](https://github.com/JuliaAstroSim/PhysicalFFT.jl) — FDM and FFT solvers for PDEs;
- [`AstroNbodySim.jl`](https://github.com/JuliaAstroSim/AstroNbodySim.jl) — N-body integration, glue layer, and parallel runtime.

By reusing these packages, developers can contribute new features without re-implementing data structures, parallel I/O, or visualization plumbing.
This makes the codebase approachable for new contributors and sustainable over the long term.

## 3. Key Features of WaveDM.jl

### 3.1 Coupled wave + N-body evolution

- A pseudo-spectral split-step Fourier solver for the 3D time-dependent SPE (§3 of the paper).
- Baryonic components can be placed **on the mesh** (`baryon_mode = :mesh`) or as **N-body particles** (`baryon_mode = :particles_static` / `:particles_dynamic`).
- Particle-based gravity supports both `DirectSum` ($\mathcal{O}(N^2)$) and `Tree` ($\mathcal{O}(N\log N)$, Gadget-2 scheme) solvers from `AstroNbodySim.jl`.
- Time-dependent tidal forces from the **Milky Way** (`MW_tidal_field = true`) and from external perturbers like the **LMC** (`LMC_tidal_field = true`) can be added.
- Baryons can be co-evolved with the wave field (bidirectional Poisson coupling).

### 3.2 Multi-level parallel execution

- **Multi-threading** — launch Julia with `julia -t N`; multi-threaded FFTW can be enabled with `FFTW.set_provider!("mkl"); FFTW.set_num_threads(Threads.nthreads())`.
- **Distributed memory** — set `distributed_memory = true` and `addprocs(N)`; WaveDM.jl uses Julia's native `Distributed.jl` (TCP/IP transport), with future MPI support planned.
- **GPU acceleration** — set `gpu = true`; WaveDM.jl dispatches to `CUDA.jl` for FFT and matrix-heavy work.

The same script runs unchanged from a laptop to a multi-node cluster.

### 3.3 Galaxy-simulation toolbox

- **Initial conditions** for gNFW / NFW / Zhao profiles, SPARC late-type galaxies, SPARC early-type galaxies, dwarf spheroidals (e.g. Crater II), and the Milky Way.
- **Trajectory lookback** of MW satellites from Battaglia et al. (2022) with optional LMC perturbations.
- **Tidal field** evaluation at the satellite's instantaneous position along its lookback orbit, with the spatial average subtracted to remove spurious center-of-mass acceleration.
- **Virial diagnostics** — total kinetic, potential, and quantum energies; radii enclosing 10 %–90 % mass fractions; velocity dispersion.
- **Real-time visualization** via `GLMakie.jl` (`Realtime = true`) and headless `UnicodePlots.jl` (`unicode_plot = true`).

## 4. Project Ecosystem

WaveDM.jl is part of the [JuliaAstroSim](https://github.com/JuliaAstroSim) ecosystem, which provides a comprehensive set of tools for astrophysical simulations:

| Package | Role |
| --- | --- |
| `AstroSimBase.jl` | Basic type definitions and unit conventions |
| `PhysicalParticles.jl` | N-body particle types and vector algebra |
| `PhysicalFFT.jl` | FFT solvers for physical problems |
| `PhysicalFDM.jl` | Finite-difference solvers for PDEs |
| `PhysicalMeshes.jl` | Mesh-based gravity solvers |
| `AstroNbodySim.jl` | N-body simulations, glue layer, parallel runtime |
| `AstroIC.jl` | Initial condition generation |
| `AstroIO.jl` | Snapshot I/O |
| `AstroPlot.jl` | Visualization and analysis |
| `GalacticDynamics.jl` | Density profiles and orbit calculations |

Together these packages form a versatile simulation stack that WaveDM.jl builds on for the SPE, N-body, and visualization layers.

## 5. Cross-Disciplinary Applications

Beyond astrophysics, the GNLSE interface and the pluggable potential $V$ let WaveDM.jl simulate a range of nonlinear Schrödinger systems:

- **Nonlinear optics** — Kerr and thermal-optical nonlinearities, beam propagation in waveguides.
- **Cold-atom physics** — Bose–Einstein condensates, vortex dynamics, dipole–dipole interactions.
- **Plasma physics** — non-relativistic wave-kinetic models.

Astrophysics, nonlinear optics, and condensed-matter physics thus share the same numerical infrastructure, enabling inter-disciplinary studies at no extra implementation cost.

## 6. When WaveDM.jl is the right tool

- ✓ You want a galaxy-scale simulation that couples a wave CDM halo to baryons, with a time-dependent tidal field.
- ✓ You want to look back the orbit of a MW satellite and run a wave-CDM simulation under a live host potential.
- ✓ You want the same code to run on a laptop, a workstation with a GPU, or a multi-node cluster.
- ✓ You want a modular framework that supports general nonlinear Schrödinger problems.

For purely cosmological wave-CDM simulations, dedicated cosmological codes (e.g. `UltraDark.jl`, `GAMER`) are more appropriate; WaveDM.jl focuses on **galactic dynamics**.
