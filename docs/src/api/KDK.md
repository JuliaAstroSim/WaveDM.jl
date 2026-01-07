# KDK

```@meta
CurrentModule = WaveDM
```

The KDK module implements the split-step Fourier method using a "kick-drift-kick" algorithm to solve the Schrödinger-Poisson equation.

## Core Functions

### apply_kick_step!

```@docs
apply_kick_step!
```

### apply_drift_step!

```@docs
apply_drift_step!
```

## Algorithm Description

The split-step Fourier method decomposes the time evolution operator into potential (kick) and kinetic (drift) components:

```math
\psi(\mathbf{r}, t + \Delta t) \approx e^{i\frac{\Delta t}{2}\hat{H}_p} e^{i\Delta t \hat{H}_k} e^{i\frac{\Delta t}{2}\hat{H}_p} \psi(\mathbf{r}, t)
```

where:
- $\hat{H}_p = \Phi$ is the potential energy operator (kick)
- $\hat{H}_k = -\frac{1}{2}\nabla^2$ is the kinetic energy operator (drift)

### Kick Step (Potential Application)

1. Compute WaveDM density from wavefunction: $\rho = |\psi|^2$
2. Calculate gravitational potential using FFT-based Poisson solver
3. Add tidal potential if enabled
4. Apply potential to wavefunction: $\psi \leftarrow e^{i\frac{\Delta t}{2}\Phi} \psi$

### Drift Step (Kinetic Evolution)

1. Transform wavefunction to Fourier space: $\tilde{\psi} = \mathcal{F}\psi$
2. Apply kinetic operator in Fourier space: $\tilde{\psi} \leftarrow e^{i\Delta t \frac{k^2}{2}} \tilde{\psi}$
3. Transform back to real space: $\psi = \mathcal{F}^{-1}\tilde{\psi}$
4. Apply boundary conditions

### Final Kick Step

Repeat the kick step to complete the second-order integration: $\psi \leftarrow e^{i\frac{\Delta t}{2}\Phi} \psi$
