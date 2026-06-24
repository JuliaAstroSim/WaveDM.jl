# KDK

```@meta
CurrentModule = WaveDM
```

The KDK module implements the split-step Fourier method using a "**kick–drift–kick**" (KDK) algorithm to evolve the Schrödinger–Poisson equation. The scheme is a second-order Strang splitting that decouples the kinetic and potential parts of the Hamiltonian.

The full mathematical derivation is given in [Algorithms §2](algorithms.md#2-the-split-step-fourier-method); this page documents the two low-level entry points that the main loop calls.

## Core Functions

### `apply_kick_step!`

```@docs
apply_kick_step!
```

The kick step **fuses steps (a), (c) and the start of (d)** of Algorithm 1 in the paper: it computes the gravitational potential from the current wavefunction, optionally adds the tidal field, applies the half-timestep potential operator, and returns the Fourier transform of the kicked wavefunction ready for the drift step.

The implementation handles:

- `device_ψ` is **consumed** by this function on the GPU path (the underlying buffer is freed before returning). The caller must not reuse the binding.
- On the GPU, intermediate buffers (`device_Φ_all`, `device_Φ_WaveDM`, `potential_grav_static`, `device_V`, `potential_all`, `nonlinear_term`, `device_ψ`) are released through `CUDA.unsafe_free!` so peak GPU memory stays close to a single wavefunction's worth.
- The tidal field (if `MW_tidal_field = true`) is added to $\Phi_{\rm total}$ at this point.

### `apply_drift_step!`

```@docs
apply_drift_step!
```

The drift step **fuses step (b) of Algorithm 1 with the absorbing boundary**: it multiplies the Fourier-space wavefunction by the linear phase factor $e^{i\,\nabla^2\,\Delta t/2}$, transforms back to real space, and applies the absorbing boundary mask $\mathcal B(\mathbf r)$.

The `spec` and `linear_phase` arguments are **consumed** by this function; their bindings should not be reused after the call.

## Algorithm Description

Write the SPE evolution operator as the sum of two pieces:

```math
\partial_t \psi \;=\; \bigl(\hat D + \hat N\bigr)\psi,
```

with the **kinetic dispersion operator**

```math
\hat D = \tfrac{i}{2}\nabla^2,
```

and the **potential operator**

```math
\hat N = -i\,\Phi_{\rm total}.
```

The second-order Strang split applied over a single timestep $\Delta t$ is

```math
\psi(\mathbf r, t+\Delta t) \;\approx\; e^{\frac{\Delta t}{2}\hat N}\,e^{\Delta t \hat D}\,e^{\frac{\Delta t}{2}\hat N}\,\psi(\mathbf r, t).
```

In Fourier space the drift is diagonal: $e^{\Delta t \hat D} \to e^{-\frac{i}{2}\Delta t\,k^2}$.

### Step (a) + (c) + start of (d) — `apply_kick_step!`

1. Compute the gravitational potential from the current density:

   ```math
   \rho = |\psi|^2, \qquad
   \Phi_{\rm DM}(\mathbf k) = -\frac{4\pi}{k^2}\bigl[\rho(\mathbf k) - \langle\rho\rangle\,\delta(\mathbf k)\bigr].
   ```

2. If `MW_tidal_field = true`, add the time-dependent MW (+LMC) tidal field to $\Phi_{\rm total}$.

3. Evaluate the user-supplied potential $V(\mathbf r, t, \psi)$ on the grid and assemble $\Phi_{\rm total}$.

4. Apply the **half-timestep potential kick** to the wavefunction:

   ```math
   \psi \;\leftarrow\; e^{\frac{\Delta t}{2}\hat N}\,\psi \;=\; e^{-i\,\Phi_{\rm total}\,\Delta t/2}\,\psi.
   ```

5. Transform to Fourier space and return the spectrum for the drift step.

### Step (b) — `apply_drift_step!`

1. Multiply the Fourier-space wavefunction by the linear phase:

   ```math
   \tilde\psi \;\leftarrow\; e^{i\,\nabla^2\,\Delta t/2}\,\tilde\psi \;=\; e^{-i\,k^2\,\Delta t/2}\,\tilde\psi.
   ```

2. Inverse Fourier transform back to real space: $\psi \leftarrow \mathcal F^{-1}\{\tilde\psi\}$.

3. Apply the **absorbing boundary mask** $\mathcal B(\mathbf r)$:

   ```math
   \psi \;\leftarrow\; \psi \odot \mathcal B(\mathbf r).
   ```

The boundary mask is set up by `setup_absorption_boundary(Xmax, Ymax, Zmax, x, y, z, absorb_coeff, dt)` and is zero (full absorption) in the outermost $\sim L/50$ cells of the box, with a smooth $\tanh$ transition.

## Sign Convention

| Operator | Definition | Evolution factor |
| --- | --- | --- |
| Kinetic $\hat D$ | $\tfrac{i}{2}\nabla^2$ | $e^{\Delta t \hat D} \to e^{-\frac{i}{2}\Delta t k^2}$ |
| Potential $\hat N$ | $-i\Phi_{\rm total}$ | $e^{\frac{\Delta t}{2}\hat N} = e^{-i\Phi_{\rm total}\Delta t/2}$ |
| Drift in code | `exp.(im .* Laplacian .* dt/2)` | matches $e^{-i k^2 \Delta t/2}$ for the convention $k \to (2\pi k_{\rm code})$ |
| Kick in code | `exp.(-im .* potential_all .* dt_kick)` | matches $e^{-i\Phi\Delta t/2}$ |

The drift implementation in [`setup_fft_operators`](@ref) uses the convention

```math
\nabla^2 \;\longrightarrow\; -(2\pi k_x)^2 - (2\pi k_y)^2 - (2\pi k_z)^2,
```

so `linear_phase = fftshift(exp.(im * Laplacian * dt/2))` evaluates to the same diagonal factor as the mathematical form above.

## Memory Management

| Path | Strategy |
| --- | --- |
| GPU (`gpu = true`) | All device buffers created in `apply_kick_step!` and `apply_drift_step!` are released through `CUDA.unsafe_free!` as soon as they are no longer needed. Peak GPU memory stays close to the size of one $\psi$ array. |
| CPU | Julia's GC reclaims host arrays when bindings go out of scope; no manual release is needed. |

The internal helper `_release!(buf, gpu)` wraps the conditional in a single function so that call sites remain readable.
