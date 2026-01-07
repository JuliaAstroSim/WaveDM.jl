# Algorithms

## 1. Schrödinger-Poisson Equation

### 1.1 Basic Formulation

The Wave Dark Matter (WaveDM) dynamics at galaxy scales are described by the time-dependent Schrödinger-Poisson equation (SPE), which governs the evolution of a macroscopic wavefunction $\psi$:

```math
\begin{align}
    i\hbar \frac{\partial\psi}{\partial t} &= -\frac{\hbar^2}{2m} \nabla^2 \psi + m\Phi \psi, \\
    \nabla^2 \Phi &= 4\pi G m |\psi|^2,
\end{align}
```

where $m$ is the boson mass, $\hbar$ is the reduced Planck constant, $G$ is the gravitational constant, and $\Phi$ is the gravitational potential. The WaveDM density is given by $\rho = m |\psi|^2$.

### 1.2 Dimensionless Form

To simplify numerical calculations, we introduce dimensionless variables following Edwards et al. (2018) and Glennon et al. (2021):

```math
\begin{aligned}
    \mathcal{L} &= \left(\frac{8\pi\hbar^{2}}{3m^{2}H_{0}^{2}\Omega_{m_{0}}}\right)^{\frac{1}{4}} \approx 121 \left(\frac{10^{-23} \mathrm{eV}}{m}\right)^{\frac{1}{2}} \mathrm{kpc}, \\
    \mathcal{T} &= \left(\frac{8\pi}{3H_{0}^{2}\Omega_{m_{0}}}\right)^{\frac{1}{2}} \approx 75.5 \mathrm{Gyr}, \\
    \mathcal{M} &= \frac{1}{G}\left(\frac{8\pi}{3H_{0}^{2}\Omega_{m_{0}}}\right)^{-\frac{1}{4}}\left(\frac{\hbar}{m}\right)^{\frac{3}{2}} \approx 7\times10^{7} \left(\frac{10^{-23} \mathrm{eV}}{m}\right)^{\frac{3}{2}} \mathrm{M}_{\odot}.
\end{aligned}
```

In these units, the SPE takes the dimensionless form (absorbing dimensional quantities for notational convenience):

```math
\begin{align}
    i \frac{\partial\psi}{\partial t} &= -\frac{1}{2} \nabla^{2} \psi + \Phi \psi, \\
    \nabla^{2} \Phi &= 4\pi |\psi|^2.
\end{align}
```

## 2. Split-step Fourier Method

### 2.1 Basic Principle

WaveDM.jl uses a pseudo-spectral split-step Fourier method to solve the SPE. This method leverages the fact that the Hamiltonian can be split into kinetic and potential parts:

```math
\hat{H} = \hat{H}_k + \hat{H}_p,
```

where $\hat{H}_k = -\frac{1}{2}\nabla^2$ is the kinetic energy operator and $\hat{H}_p = \Phi$ is the potential energy operator. The time evolution operator can be approximated using the Baker-Campbell-Hausdorff (BCH) formula:

```math
\psi(\mathbf{r}, t+\Delta t) \approx e^{i\Delta t \hat{H}} \psi(\mathbf{r}, t) \approx e^{i\frac{\Delta t}{2}\hat{H}_p} e^{i\Delta t \hat{H}_k} e^{i\frac{\Delta t}{2}\hat{H}_p} \psi(\mathbf{r}, t).
```

### 2.2 Computational Steps

The split-step method proceeds in three stages:

1. **Potential Kick (First Half)**:
   Apply the potential operator in real space:
   ```math
   \psi_1 = e^{i\frac{\Delta t}{2}\hat{H}_p} \psi(\mathbf{r}, t) = e^{i\frac{\Delta t}{2}\Phi} \psi(\mathbf{r}, t).
   ```

2. **Kinetic Drift**:
   Apply the kinetic operator in Fourier space:
   ```math
   \psi_2 = \mathcal{F}^{-1} \left\{ e^{i\Delta t \hat{H}_k} \mathcal{F} \{ \psi_1 \} \right\},
   ```
   where $\mathcal{F}$ denotes the Fourier transform and $e^{i\Delta t \hat{H}_k}$ is diagonal in Fourier space:
   ```math
   e^{i\Delta t \hat{H}_k} = e^{i\frac{\Delta t}{2} k^2}, 
   ```
   with $\mathbf{k}$ being the wavevector.

3. **Potential Kick (Second Half)**:
   Apply the potential operator again in real space:
   ```math
   \psi(\mathbf{r}, t+\Delta t) = e^{i\frac{\Delta t}{2}\hat{H}_p} \psi_2 = e^{i\frac{\Delta t}{2}\Phi} \psi_2.
   ```

### 2.3 Poisson Equation Solving

The gravitational potential $\Phi$ is obtained by solving the Poisson equation:

```math
\nabla^2 \Phi = 4\pi \rho.
```

WaveDM.jl uses a Fast Fourier Transform (FFT) based Poisson solver, which transforms the equation to Fourier space:

```math
-k^2 \Phi(\mathbf{k}) = 4\pi \rho(\mathbf{k}),
```

resulting in:

```math
\Phi(\mathbf{k}) = -\frac{4\pi}{k^2} \rho(\mathbf{k}).
```

To avoid the singularity at $k=0$, we subtract the mean density:

```math
\Phi(\mathbf{k}) = -\frac{4\pi}{k^2} (\rho(\mathbf{k}) - \langle\rho\rangle \delta(\mathbf{k})),
```

where $\langle\rho\rangle$ is the spatial average of $\rho$.

## 3. Numerical Implementation

### 3.1 Spatial Discretization

WaveDM.jl uses a uniform Cartesian grid with $N_x \times N_y \times N_z$ grid points. The coordinates of a grid point $(i,j,k)$ are:

```math
\begin{aligned}
    x_i &= -\frac{L_x}{2} + i \Delta x, \\
    y_j &= -\frac{L_y}{2} + j \Delta y, \\
    z_k &= -\frac{L_z}{2} + k \Delta z,
\end{aligned}
```

where $L_x, L_y, L_z$ are the box dimensions and $\Delta x, \Delta y, \Delta z$ are the grid spacings.

### 3.2 Boundary Conditions

WaveDM.jl supports two types of boundary conditions:

1. **Periodic Boundaries**:
   ```math
   \psi(x+L_x, y, z, t) = \psi(x, y, z, t), 
   ```
   Suitable for cosmological simulations or isolated systems with large boxes.

2. **Absorbing Boundaries**:
   A damping term is added at the boundaries to absorb outgoing waves:
   ```math
   \psi(x, y, z, t+\Delta t) = e^{-\gamma(x,y,z)\Delta t} \psi(x, y, z, t),
   ```
   where $\gamma(x,y,z)$ is a damping function that increases near the boundaries.

### 3.3 Time Integration

The time step $\Delta t$ is determined by the Courant-Friedrichs-Lewy (CFL) condition:

```math
\Delta t \leq 0.5 \min\left( \frac{4}{3\pi} \Delta x^2, \frac{2\pi}{\max|\Phi|} \right),
```

where the first term ensures stability of the kinetic term and the second term limits phase variations from the potential.

### 3.4 GPU Acceleration

WaveDM.jl leverages CUDA for GPU acceleration. Key operations that benefit from GPU computing include:

- Fast Fourier Transforms (FFT)
- Element-wise operations on large arrays
- Gradient calculations
- Potential calculations

The code automatically detects available GPUs and uses them if enabled by the user.

## 4. Initial Conditions

### 4.1 Density Profiles

WaveDM.jl supports several dark matter density profiles for initial conditions:

1. **Generalized Navarro-Frenk-White (gNFW)**:
   ```math
   \rho(r) = \frac{\rho_0}{(r/r_s)^\beta (1 + r/r_s)^{3-\beta}}.
   ```

2. **Zhao Profile**:
   ```math
   \rho(r) = \frac{\rho_0}{(r/r_s)^\gamma (1 + (r/r_s)^\alpha)^{(3-\gamma)/\alpha}}.
   ```

3. **Navarro-Frenk-White (NFW)**:
   ```math
   \rho(r) = \frac{\rho_0}{(r/r_s) (1 + r/r_s)^2}.
   ```

### 4.2 Velocity Field Initialization

The initial velocity field is constructed using a combination of circular rotation and random components:

```math
\mathbf{v} = \zeta_{\mathrm{rot}} |v_c| \hat{\mathbf{e}}_\phi + (1-\zeta_{\mathrm{rot}}) |v_c| \hat{\mathbf{e}}_{\mathrm{random}},
```

where $v_c = \sqrt{r\nabla\Phi}$ is the circular velocity, $\hat{\mathbf{e}}_\phi$ is the azimuthal unit vector, and $\zeta_{\mathrm{rot}}$ controls the fraction of ordered rotation.

## 5. Gravitational Potential Solver

### 5.1 WaveDM Component

The WaveDM contribution to the potential is computed using the FFT-based Poisson solver described in Section 2.3.

### 5.2 Baryonic Component

For baryonic components, WaveDM.jl uses N-body methods to compute the gravitational potential:

1. **Direct Summation (DS)**:
   ```math
   \Phi_b(\mathbf{r}) = -G \sum_{i=1}^{N} \frac{m_i}{|\mathbf{r}-\mathbf{r}_i|},
   ```
   Suitable for moderate particle counts.

2. **Tree-based Solver**:
   Uses an octree structure to compute forces with $O(N\log N)$ complexity, following the Gadget-2 algorithm.

## 6. Tidal Field Integration

WaveDM.jl includes time-dependent tidal fields to simulate the evolution of Milky Way satellites:

1. **Milky Way Tidal Field**:
   Computed from high-resolution trajectory lookback data.

2. **Large Magellanic Cloud (LMC) Influence**:
   Can be included as a time-varying perturber with its own trajectory.

## 7. Visualization and Analysis

### 7.1 Real-time Visualization

WaveDM.jl provides real-time visualization of:
- Density fields
- Gravitational potentials
- Velocity fields
- Rotation curves

### 7.2 Post-processing

After simulation, WaveDM.jl generates:
- Density profiles
- Rotation curves
- Power spectra
- Lagrangian radii evolution
- Virial energy evolution

### 7.3 Profile Fitting

The code automatically fits density profiles to the simulated halos, allowing comparison with observational data.
