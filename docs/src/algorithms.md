# Algorithms

This page describes the numerical methods used by WaveDM.jl. It mirrors the structure of §§2–3 of the accompanying paper and uses the same sign conventions. All equations are written in the dimensionless units introduced in [Introduction §1.3](introduction.md#13-dimensionless-units).

## 1. The Schrödinger–Poisson System

### 1.1 Dimensionless form

```math
\begin{aligned}
i\,\partial_t \psi(\mathbf r, t) &= -\tfrac{1}{2}\nabla^2\psi + \Phi_{\rm total}(\mathbf r, t)\,\psi, \\
\nabla^2 \Phi_{\rm total}(\mathbf r, t) &= 4\pi\bigl(\rho_{\rm DM} + \rho_{\rm baryon}\bigr),
\end{aligned}
```

with $\rho_{\rm DM} = |\psi|^2$ and $\Phi_{\rm total} = \Phi_{\rm DM} + \Phi_{\rm baryon}$. The wavefunction is written in polar form as $\psi = \sqrt{\rho}\,e^{i\theta}$, with velocity field $\mathbf v = \nabla\theta$.

### 1.2 Generalized form

For non-astrophysical applications WaveDM.jl solves the **generalized nonlinear Schrödinger equation (GNLSE)**:

```math
i\,\partial_t \psi(\mathbf r, t) \;=\; -\tfrac{1}{2}\nabla^2\psi \;+\; V(\mathbf r, t, \psi)\,\psi,
```

where $V$ is an arbitrary, user-supplied potential. The two native classes are:

- **Local** — $V = F(|\psi(\mathbf r, t)|)$, depending only on the field amplitude at the same point. Example: Kerr nonlinearity $V = \kappa|\psi|^2$ used for self-interacting wave CDM and instantaneous nonlinear optical response.
- **Nonlocal** — $V$ is an integral over the field distribution, $V = \int U(\mathbf r - \mathbf r')\,|\psi(\mathbf r')|^2 d^3 r'$, or equivalently the solution of an auxiliary field equation (e.g. Poisson for self-gravity, dipole–dipole in atomic condensates, thermal-optical nonlinearities in NLO).

### 1.3 Madelung form and the quantum energy

The Madelung transformation yields the continuity equation and the quantum Bernoulli equation. The *quantum potential*

```math
Q = -\frac{\hbar^2}{2m}\,\frac{\nabla^2\sqrt\rho}{\sqrt\rho}
   = -\frac{\hbar^2}{4m}\left[\frac{\nabla^2\rho}{\rho} - \frac{1}{2}\frac{(\nabla\rho)^2}{\rho^2}\right]
```

encodes the Heisenberg uncertainty principle. Its density-weighted volume integral is the **quantum (gradient) energy** $\mathcal Q$:

```math
\mathcal Q \;\equiv\; \int \rho Q\,d^3 r
             \;=\; \frac{\hbar^2}{2m}\!\int\! |\nabla\sqrt\rho|^2 d^3r,
```

where the surface term vanishes for any realistic NFW-like halo (power-law decay $\rho \propto r^{-n}$ with $n > 1$). In dimensionless units $\mathcal Q = \int |\nabla\sqrt\rho|^2 d^3r$.

The **quantum pressure** is $-(1/m)\nabla Q$ and appears in the quantum Euler equation. Three synonyms appear in the literature: *quantum energy* (Hui et al. 2017), *gradient energy* $K_\rho$ (Mocz et al. 2017), and *effective thermal energy* (Liao et al. 2025) — they are all the same density-weighted volume integral of $Q$.

The total energy of an isolated, fully self-gravitating halo is

```math
\mathcal E = \mathcal K + \mathcal V + \mathcal Q,
```

where $\mathcal K$ is the volume-integrated kinetic energy and $\mathcal V$ the gravitational potential energy. In simulations with fixed baryonic backgrounds and absorbing boundaries $\mathcal E$ is not strictly conserved (see Section 5.1 below).

For a wave-CDM halo governed by the SPE, the **virial theorem** at virial equilibrium takes the form

```math
2\,\mathcal K + 2\,\mathcal Q + \mathcal V = 0,
```

which can be evaluated from a single simulation snapshot (or, equivalently, from single-epoch observational data of a realistic galaxy). This virial formula is the most stringent test available without time-series information.

## 2. The Split-Step Fourier Method

### 2.1 Operator splitting

Write the evolution operator as the sum of two pieces:

```math
\partial_t \psi \;=\; \bigl(\hat D + \hat N\bigr)\psi,
```

with the **kinetic dispersion operator** $\hat D = \tfrac{i}{2}\nabla^2$ and the **potential operator** $\hat N = -i\,\Phi_{\rm total}$. Both pieces are exactly solvable (linear in their respective representation). The Strang second-order symmetric split is

```math
\psi(\mathbf r, t+\epsilon) \;\approx\; e^{\frac\epsilon 2 \hat N}\,e^{\epsilon \hat D}\,e^{\frac\epsilon 2 \hat N}\,\psi(\mathbf r, t).
```

This is the "**kick–drift–kick**" (KDK) time-integration scheme. In Fourier space the drift is diagonal:

```math
e^{\epsilon \hat D}\;\longrightarrow\;e^{-\frac{i}{2}\epsilon\,k^2}.
```

### 2.2 The four computational steps

Let $\mathcal F$ denote the forward FFT and $\mathcal F^{-1}$ its inverse. With the time-averaging potential update between the two kicks, the KDK step becomes:

```math
\begin{aligned}
&\text{(a)} &\psi_1 &\;\leftarrow\; e^{\frac\epsilon 2\hat N}\psi(\mathbf r, t), \\
&\text{(b)} &\tilde\psi &\;\leftarrow\; \mathcal F\{\psi_1\},\quad
                        \tilde\psi \;\leftarrow\; e^{-\frac{i}{2}\epsilon k^2}\,\tilde\psi,\quad
                        \psi_2 \;\leftarrow\; \mathcal F^{-1}\{\tilde\psi\}, \\
&\text{(c)} &\rho &\;\leftarrow\; |\psi_2|^2, \\
&          &\Phi_{\rm DM} &\;\leftarrow\; \mathcal F^{-1}\!\left\{-\frac{4\pi}{k^2}\mathcal F\bigl\{\rho - \langle\rho\rangle\bigr\}\right\}, \\
&          &\Phi_{\rm total} &\;\leftarrow\; \Phi_{\rm DM} + \Phi_{\rm baryon}, \\
&\text{(d)} &\psi(\mathbf r, t+\epsilon) &\;\leftarrow\; e^{\frac\epsilon 2\hat N}\,\psi_2, \\
&          &\text{absorbing boundary:} &\quad \psi \;\leftarrow\; \psi \odot \mathcal B(\mathbf r), \\
&          &\text{snapshot if }n &\bmod n_{\rm stride} = 0.
\end{aligned}
```

Step (c) — the gravity update — is what makes the scheme a *self-consistent* SPE solver: the kinetic drift uses the old potential, but the second kick uses an updated potential that incorporates the new density field. This averaging is essential for second-order accuracy.

### 2.3 Algorithm 1 — split-step Fourier loop

The full main loop, in pseudocode, is:

```text
input : ψ(r, 0), Δt, T_max
output: ψ(r, T_max)
n_max      ← ⌊T_max / Δt⌋
n_stride   ← ⌊T_snap / Δt⌋

for n = 0, 1, …, n_max − 1 do
    # Step (a): first half potential kick
    ψ ← exp(Δt/2 · N̂) ψ        with N̂ = −i Φ_total

    # Step (b): kinetic drift in Fourier space
    ψ̃ ← F{ψ}                   # forward FFT
    ψ̃ ← exp(−i/2 · Δt · k²) · ψ̃
    ψ ← F⁻¹{ψ̃}                 # inverse FFT

    # Step (c): update gravitational potential
    ρ  ← |ψ|²
    Φ_DM     ← F⁻¹{ −4π/k² · F{ρ − ⟨ρ⟩} }
    Φ_total  ← Φ_DM + Φ_baryon

    # Step (d): second half potential kick
    ψ ← exp(Δt/2 · N̂) ψ

    # Absorbing boundary condition
    ψ ← ψ ⊙ B(r)

    if n mod n_stride = 0 then save ψ(r, t = n·Δt)
end
return ψ(r, T_max)
```

In code the steps (a), (b), (c), (d) map to [`apply_kick_step!`](@ref) and [`apply_drift_step!`](@ref) in the KDK module — see [KDK API](@ref) for signatures.

### 2.4 FFT complexity

$\mathcal F$ and $\mathcal F^{-1}$ use `FFTW.jl` (CPU) or `CUDA.jl` (GPU). Both run in $\mathcal O(N\log N)$ for a grid with $N$ cells, so an entire KDK step is $\mathcal O(N\log N)$.

## 3. Spatial Discretization

### 3.1 Uniform Cartesian grid

The simulation uses a uniform Cartesian grid with $N_x\times N_y\times N_z$ points (keywords `Nx`, `Ny`, `Nz`). For nearly spherical halos a cubic grid $N^3$ is conventional. The cell sizes are

```math
\Delta x = \frac{L_x}{N_x-1}, \quad
\Delta y = \frac{L_y}{N_y-1}, \quad
\Delta z = \frac{L_z}{N_z-1},
```

and grid coordinates are

```math
x_i = -\tfrac{L}{2} + i\,\Delta x,\; i = 0,\dots,N_x-1,
```

with the same form for $y_j$, $z_k$ (keywords `Xmax`, `Ymax`, `Zmax` set the half-side lengths; the full side length is $L = 2X_{\max}$). Gradient operators on the grid are computed with a second-order central difference stencil.

### 3.2 Resolution criteria

Five criteria guide the choice of grid resolution. The code prioritises criterion (1) and dynamically adjusts $L_x$ accordingly.

| # | Criterion | Reason |
| --- | --- | --- |
| 1 | $\Delta x \lesssim \lambda_{\rm dB}/3$ | Resolves quantum oscillations and avoids phase aliasing |
| 2 | $\Delta x \lesssim r_s/3$ | Adequately samples the density field around the scale radius |
| 3 | $\Delta x \lesssim 0.3$ kpc | Accurate inner-slope fitting for dwarf halos |
| 4 | $L_x \gtrsim 10\,\lambda_{\rm dB}$ | Suppresses boundary artifacts in the soliton core |
| 5 | $L_x \gtrsim 8\,r_s$ | Adequate gravitational potential resolution |

For $m_{22} \equiv m_a / 10^{-22}\,{\rm eV} = 1.0$ and $v_{\rm rot}\sim 130$ km s$^{-1}$ one obtains $\lambda_{\rm dB} \sim 0.77$ kpc, giving the practical recommendations $\Delta x \le 0.4$ kpc and $N \ge 256$ for $L = 80$ kpc. Empirically, $\Delta x \lesssim \lambda_{\rm dB}/3$ is needed for unbiased velocity fields, requiring $N \ge 512$ at $m_{22} = 1.0$. Grid dimensions should be powers of two for FFT efficiency.

The maximum resolvable velocity is $v_{\max} = (\hbar/m)(\pi/\Delta x)$; the conservative choice $v_{\max} \ge 2\max\{v_{\rm rot}\}$ is required for unbiased dynamics.

## 4. Boundary Conditions

### 4.1 Periodic (default)

`boundary = Periodic()` is the natural choice for the Fourier-based solver. The Poisson equation is solved in Fourier space by

```math
\Phi_{\rm DM}(\mathbf k) = -\frac{4\pi}{k^2}\bigl[\rho(\mathbf k) - \langle\rho\rangle\,\delta(\mathbf k)\bigr],
```

where $\langle\rho\rangle$ is the spatial average of the density — subtracting it removes the $\mathbf k = 0$ singularity of the Laplacian.

To approximate isolated (vacuum) boundary conditions on a periodic domain, the box size should be at least eight times the halo scale radius, ensuring that no substantial bulk mass crosses the boundaries. This is the recommended choice for production runs.

### 4.2 Vacuum (isolated)

`boundary = Vacuum()` uses the **four-step algorithm of James (1977)** for isolated systems. It is $\mathcal O(N^4)$ for an $N^3$ grid, so it is best reserved for **verification** of periodic-boundary runs at modest resolution.

### 4.3 Absorbing layer

For simulations where wave energy must exit the domain without reflection, WaveDM.jl multiplies $\psi$ by a smooth absorbing mask $\mathcal B(\mathbf r)$ at the end of every KDK step:

```math
\mathcal B(\mathbf r) = \exp\!\left[-\alpha\!\left(
    6
    - \tanh\!\frac{x_i + L_x/2}{L_x/N_b}
    + \tanh\!\frac{x_i - L_x/2}{L_x/N_b}
    - \tanh\!\frac{y_j + L_y/2}{L_y/N_b}
    + \tanh\!\frac{y_j - L_y/2}{L_y/N_b}
    - \tanh\!\frac{z_k + L_z/2}{L_z/N_b}
    + \tanh\!\frac{z_k - L_z/2}{L_z/N_b}
\right)\right]\Delta t.
```

The defaults are $\alpha = 10$ (keyword `absorb_coeff`) and $N_b = 50$ (transition scale $\sim L/50 \sim \mathcal O(10^2)$ cells). The mask limits unphysical interference patterns within the simulated halo.

## 5. Time Integration

### 5.1 Total energy drift

WaveDM.jl tracks three volume-integrated energies:

- **Kinetic** — $\mathcal K = \frac{1}{2}\int (\rho\mathbf v)^2 / \rho\,d^3 r$, where the momentum density is $\rho\mathbf v = (\psi^*\nabla\psi - \psi\nabla\psi^*)/(2i)$.
- **Gravitational potential** — $\mathcal V = \int \rho\,\Phi\,d^3 r$.
- **Quantum gradient** — $\mathcal Q = \int |\nabla\sqrt\rho|^2 d^3 r$.

The total energy $\mathcal E = \mathcal K + \mathcal V + \mathcal Q$ is *approximately* conserved for a virialized, isolated halo. In a fixed baryonic background with absorbing boundaries, however, $\mathcal E$ drifts downward (mass loss through the absorber) — this is a real physical effect, not a numerical artefact, and is monitored in real time.

### 5.2 Timestep criterion

Following Schwabe et al. (2016), Mocz et al. (2017), May & Springel (2021) and May et al. (2023), WaveDM.jl uses

```math
\Delta t \;\le\; \Delta t_{\max} \;=\; 0.5\,\min\!\left(
    \frac{4}{3\pi}\frac{m}{\hbar}\Delta x^2,\;
    2\pi\frac{\hbar}{m}\frac{1}{\max|\Phi_c|},\;
    2\pi\frac{\hbar}{m}\frac{1}{|\kappa|\max\{|\psi|^2\}}
\right),
```

where the three terms limit phase variations from the **kinetic propagation**, the **gravitational potential**, and the **self-interaction** respectively. In dimensionless units the prefactors $m/\hbar$ are absorbed and the expression reads

```math
\Delta t_{\max} = 0.5\,\min\!\left(
    \frac{4}{3\pi}\Delta x^2,\;
    \frac{2\pi}{\max|\Phi_c|},\;
    \frac{2\pi}{|\kappa|\max|\psi|^2}
\right).
```

Pass `autoset_timestep = true` to use it; the safety factor $\zeta_t$ is set by `autoset_timestep_ratio` (default 0.9).

### 5.3 Simulation duration

The ground-state oscillation period is $\tau_{00} \sim \hbar/(m_a v_{\rm vir}^2)$ [Chiang et al. 2021]. For Crater II ($v_{\rm vir}\approx 10$ km s$^{-1}$, $m_{22}=50.0$) one obtains $\tau_{00} \approx 3.75$ Myr. Simulations are typically run for several Gyr so that virial equilibrium is reached (≈ 1.5 Gyr for Crater II). Setting `Tmax = 5u"Gyr"` is a safe default for halos of Crater-II mass.

## 6. Quasi-Stationary Equilibrium (QSS) Diagnostics

Three nested criteria are provided, in order from global to spatially resolved.

### 6.1 Overall energy terms (global, cheap)

Monitor $\mathcal K(t)$, $\mathcal V(t)$, $\mathcal Q(t)$ and the second time derivative of the moment of inertia. After virialization, $\mathcal E$ approaches a constant and $d^2 I / dt^2 \to 0$.

**Use when:** the analysis only needs integrated halo properties; you want the cheapest indicator.

### 6.2 Mass-fraction radii (spatial, cheap)

Track the radii enclosing fixed mass fractions (10 %, 20 %, …, 90 %). These radii initially increase as mass redistributes outward during relaxation, then approach steady values. Stabilization of all nine curves indicates that the spatial configuration of $\rho(r)$ no longer evolves.

**Use when:** the analysis requires that $\rho(r)$ be approximately stationary in time, e.g. when fitting a single density profile.

### 6.3 Velocity-dispersion profile (post-processing, expensive)

For a halo in a quasi-stationary state supported by both velocity dispersion and rotational motion, the Jeans equation

```math
\frac{d(\rho\sigma_r^2)}{dr} + \frac{2\beta(\rho\sigma_r^2)}{r} + \rho\frac{\langle v_\theta\rangle^2 + \langle v_\phi\rangle^2}{r} = -\rho\,\frac{d\Phi_{\rm total}}{dr} - \rho\,\frac{dQ}{dr},
```

where $\beta = 1 - (\sigma_\theta^2 + \sigma_\phi^2)/(2\sigma_r^2)$ is the velocity-anisotropy parameter and $Q$ is the quantum potential. The right-hand side represents the gravitational force and the quantum-pressure force; the three terms on the left correspond to the radial pressure-gradient term, the anisotropic-pressure correction, and centrifugal support from rotational streaming.

Assuming isotropic velocity dispersion ($\beta \equiv 0$) and rearranging gives the predicted radial velocity dispersion

```math
\sigma_r^2(r) = \frac{1}{\rho(r)} \int_r^\infty \rho(r')\left(\frac{d\Phi_{\rm total}}{dr'} + \frac{dQ}{dr'} + \frac{\langle v_\theta\rangle^2 + \langle v_\phi\rangle^2}{r'}\right) dr'.
```

For an isolated halo this Jeans-derived $\sigma_r(r)$ matches the directly measured radial velocity dispersion in all three spherical components, validating the QSS. With a time-dependent tidal field the steady-state assumption breaks down; the simulated and Jeans-predicted $\sigma_r$ no longer agree, signalling a *tide-driven non-equilibrium state*.

**Use when:** the science goal is to compare with integral-field-spectroscopy observations; this is a post-processing diagnostic because it requires the velocity moments on the grid.

## 7. Initial Conditions

### 7.1 Density profile

The macroscopic amplitude derives from a target density profile:

```math
|\psi| = \sqrt{\rho_{\rm halo} / m},
```

with the density profile chosen from `gNFW`, `NFW`, `Zhao`, or a user-supplied function:

- **gNFW** — $\rho(r) = \rho_0 \, / \bigl[(r/r_s)^\beta (1 + r/r_s)^{3-\beta}\bigr]$.
- **NFW** — $\beta = 1$ special case of gNFW.
- **Zhao** — $\rho(r) = \rho_0 \, / \bigl[(r/r_s)^\gamma (1 + (r/r_s)^\alpha)^{(3-\gamma)/\alpha}\bigr]$.

### 7.2 Velocity field

The phase $\theta(\mathbf r)$ encodes the initial velocity field through $\mathbf v = \nabla\theta$. WaveDM.jl uses **direct integration** of the velocity components by default. An alternative is to solve the Poisson equation $\nabla^2\theta = \nabla\cdot\mathbf v$ with homogeneous Dirichlet boundaries, which guarantees a curl-free field.

The velocity magnitude is taken from the local circular velocity $v_c = \sqrt{r\,d\Phi_{\rm total}/dr}$ and assembled as a weighted sum of rotation and random motion:

```math
\mathbf v = \zeta_{\rm vel}\bigl[\zeta_{\rm rot}\,v_c\,\hat{\mathbf e}_\phi
                                + (1-\zeta_{\rm rot})\,v_c\,\hat{\mathbf e}_{\rm random}\bigr],
```

with $\zeta_{\rm vel}$ controlled by `velocity_ratio` and $\zeta_{\rm rot}$ by `rotational_ratio`.

### 7.3 Pooling strategy (velocity coarse-graining)

A fully random velocity field on each cell typically has too small a divergence, leading to underestimated phase gradients, repeated halo collapses, and unrealistic density distributions. To prevent this, WaveDM.jl implements a **pooling** strategy (`bulk_perturb = true`, default `bulk_size = 4`): within each $4\times 4\times 4$ block of cells the velocity vectors are reset to the value of a reference cell — either the block's corner, the box centre, or a randomly placed cell. The pooling preserves velocity coherence without introducing spurious gradients.

For purely rotational velocity fields ($\zeta_{\rm rot}=1$) pooling is unnecessary; for mixed configurations it is a small but important correction.

### 7.4 Baryonic initial conditions

Baryons are placed either on the same mesh (`baryon_mode = :mesh`) or as N-body particles (`baryon_mode = :particles_static` / `:particles_dynamic`). The Milky Way mass model used in the case study follows Zhu et al. (2023):

```math
\begin{aligned}
\rho_{\rm bulge}(r') &= \frac{\rho_{0,b}}{(1+r'/r_0)^\alpha}\exp\!\bigl[-(r'/r_{\rm cut})^2\bigr], \\
\rho_{\rm disk}(R, z) &= \frac{\Sigma_0}{2 z_d}\exp\!\bigl(-|z|/z_d - R/R_d\bigr), \\
\rho_{\rm gas}(R, z) &= \frac{\Sigma_0}{4 z_d}\exp\!\bigl(-R_m/R - R/R_d\bigr)\,{\rm sech}^2(z/2z_d), \\
\rho_{\rm DM}(r) &= \rho_{0,h}\!\left(\frac{r}{r_h}\right)^{-\gamma}
                 \!\left[1 + \left(\frac{r}{r_h}\right)^\alpha\right]^{(\gamma-\beta)/\alpha},
\end{aligned}
```

with $r' = \sqrt{R^2 + (z/q)^2}$ in cylindrical coordinates. Sampling the baryons as N-body particles is illustrated in [Examples](@ref).

## 8. Tidal Field

When `MW_tidal_field = true`, the gravitational potential of the MW (and, optionally, the LMC with `LMC_tidal_field = true`) is added to $\Phi_{\rm total}$ at every KDK step. The tidal field is evaluated at the satellite's instantaneous position $\mathbf r_{\rm sat}(t)$ along its lookback orbit (Section 9), so the wave-CDM halo experiences a time-dependent external potential.

To prevent the tidal field from exerting a spurious net force on the halo, WaveDM.jl subtracts the spatial average of the tidal gradient over the computational domain (equivalently, the mean acceleration). Without this correction the halo would acquire unphysical center-of-mass motion.

## 9. Trajectory Lookback

Trajectory lookback uses the positional and velocity data of MW satellites from Battaglia et al. (2022). The tool converts (RA, Dec, $\mu_{\alpha,*}$, $\mu_\delta$, distance, radial velocity) to the Galactocentric frame and integrates the orbit backward in the (possibly time-dependent) MW potential.

LMC perturbations are included by integrating the LMC orbit backward and adding its gravitational influence to the satellite's equation of motion. The LMC is modeled with $M_{\rm LMC} = 1.5\times 10^{11}\,M_\odot$ and $r_s = 10.84$ kpc, and can be either an analytic Plummer potential or a system of N-body particles (`particles_LMC`).

The resulting trajectory $\mathbf r_{\rm sat}(t)$ drives the time-dependent tidal field in Section 8. See the [Examples](@ref) page for a Crater II script.

## 10. Convergence Test Strategy

The following four axes are exercised in the validation suite (see [Tests](https://github.com/JuliaAstroSim/WaveDM.jl/tree/main/test)):

- **Timestep $\Delta t$** — vary by factors of 0.5× and 2×. Should be negligible impact on density and velocity-dispersion profiles.
- **Spatial resolution $\Delta x$** — vary by factors of 0.75× and 1.5×. Coarser grids lock the density into its initial distribution; finer grids improve interference-pattern resolution.
- **Domain size $L$** — vary by factors of 0.67× and 1.33×. Boundary artifacts become visible at small $L/r_s$.
- **Initial condition parameters** — $\zeta_{\rm rot}$ (rotation), $\zeta_{\rm pool}$ (`bulk_size`), and the absorber coefficient $\alpha$.

The reference simulation uses $m_{22}=50$, $N_x=384$, $L=40$ kpc, $N_t=377$, $\zeta_{\rm rot}=0.9$, $\alpha=10$, $\zeta_{\rm pool}=4$, and $T_{\max}=3$ Gyr. Tests are run with the MW tidal field **disabled** so the halo is free to virialize cleanly.
