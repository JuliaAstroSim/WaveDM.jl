# Examples

This page walks through the end-to-end simulations discussed in §3 of the accompanying paper. Each example is a complete, runnable Julia script. Run from the Julia REPL or `julia --project` with the indicated packages installed.

## 0. Smoke test (CPU, 64³)

The fastest way to verify that WaveDM.jl is correctly installed and that the SPE loop is wired up. Runs in a few seconds on any machine.

```julia
using WaveDM

simulate_waveDM(;
    model = :MW,
    Xmax = 2.0,
    Tmax = 0.001,
    Nt   = 3,
    Nx = Ny = Nz = 64,
    gpu = false,
    Realtime = false,
    unicode_plot = true,             # terminal progress
    outputdir = mktempdir(),
)
```

You should see several `@info` lines, a terminal heat-map of the density slice, and no errors.

## 1. Isolated wave-DM halo (free propagation, Gaussian IC)

Useful for verifying the KDK integrator in the simplest possible setup — no self-gravity, no baryons.

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model = :MW,
    IC    = (x, y, z) -> exp(-(x^2 + y^2 + z^2) / 0.1^2),
    V     = (x, y, z, ψ) -> 0.0,
    κ     = 0,
    Nx = Ny = Nz = 128,
    Xmax = 2.0,
    Tmax = 0.01,
    Nt   = 200,
    autoset_timestep = true,
    gpu = false,
    Realtime = true,
    title = "FreeProp",
    outputdir = "results/FreeProp",
)
```

The wave packet disperses under the kinetic operator; with `autoset_timestep = true` the timestep is chosen from the CFL criterion.

## 2. Milky Way ultra-faint dwarf — Crater II

The flagship case study from §3 of the paper. Crater II is modelled with a gNFW density profile:

```math
\rho_{\rm DM}(r) = \frac{\rho_0}{(r/r_s)^\beta (1 + r/r_s)^{3-\beta}},
```

with the Crater-II parameters from Hayashi et al. (2023). The boson mass is $m_{22} = 50.0$, the resolution is $384^3$, the box half-side is 20 kpc, and the run is 6 Gyr long. Baryons are neglected (Crater II is DM-dominated).

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model       = :dwarf_UFDs,
    Galaxy_id   = 6,                          # Crater II
    V           = (x, y, z, ψ) -> 0.0,
    Nx          = 384,
    Xmax        = 20u"kpc",
    Tmax        = 6.0u"Gyr",
    autoset_timestep = true,
    absorb_coeff     = 10,                    # moderate absorbing layer
    gpu         = true,
    Realtime    = true,
    plot_virial = true,                       # monitor virial energies
    title       = "CraterII",
    outputdir   = "results/CraterII",
)
```

The simulation writes an overview figure, an MP4 of the live density slice, and a CSV with the time series of total mass, mass-fraction radii, and virial energies.

### 2.1 Crater II with the Milky Way tidal field

To add a time-dependent MW tidal field (and, optionally, the LMC perturbation) to the same halo:

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model       = :dwarf_UFDs,
    Galaxy_id   = 6,
    V           = (x, y, z, ψ) -> 0.0,
    Nx          = 384,
    Xmax        = 20u"kpc",
    Tmax        = 6.0u"Gyr",
    autoset_timestep = true,
    absorb_coeff     = 10,
    MW_tidal_field = true,                    # MW tidal field
    LMC_tidal_field = true,                   # LMC perturbation
    gpu         = true,
    title       = "CraterII_MW_LMC",
    outputdir   = "results/CraterII_MW_LMC",
)
```

The halo ends up in a **tide-driven non-equilibrium state** — the radially averaged velocity dispersion no longer matches the Jeans prediction, which is a useful diagnostic for stripping studies (see [Algorithms §6.3](algorithms.md#63-velocity-dispersion-profile-post-processing-expensive)).

## 3. Sampling the Milky Way as a static N-body background

If you want to drive a wave-CDM simulation with an explicit N-body representation of the MW, sample the baryons with `AstroIC` and pass them in via `baryon_mode = :particles_static`:

```julia
using WaveDM
using AstroIC
using PhysicalParticles, UnitfulAstro, Unitful

# Thin stellar disk
disk_thin = generate(
    AstroIC.ExponentialDisc(;
        collection       = STAR,
        NumSamples       = 50_000,
        TotalMass        = 3.5e10u"Msun",
        ScaleRadius      = 2.42u"kpc",
        ScaleHeight      = 0.30u"kpc",
    )
)

# HI gas disk
gas_HI = generate(
    AstroIC.ExponentialDisc(;
        collection       = GAS,
        NumSamples       = 30_000,
        TotalMass        = 8.0e9u"Msun",
        ScaleRadius      = 7.0u"kpc",
        ScaleHeight      = 0.085u"kpc",
        HoleRadius       = 4.0u"kpc",         # R_m in the Zhu 2023 model
    )
)

# Add a bulge, thick disk, etc. — see the paper for the full MW model
particles = [disk_thin; gas_HI]

simulate_waveDM(;
    model              = :MW,
    baryon_mode        = :particles_static,
    baryon_particles   = particles,
    GravitySolver      = Tree(),
    SofteningLength    = 1.0u"kpc",
    Nx = Ny = Nz = 256,
    Xmax = 100u"kpc",
    Tmax = 3.0u"Gyr",
    autoset_timestep = true,
    gpu = true,
    title = "MW_static_Nbody",
    outputdir = "results/MW_static_Nbody",
)
```

`baryon_mode = :particles_static` means the MW particles are held fixed — their gravitational potential is evaluated at every SPE step but they do not move. This is the cheapest way to add a realistic host potential.

## 4. Trajectory lookback of Crater II

`AstroIC.load_data_MW_satellites()` returns the Battaglia et al. (2022) observational catalogue. The script below integrates the Crater II orbit backward 6 Gyr in a time-dependent MW potential, with the LMC perturbation included. The resulting trajectory can then be fed to the wave-DM simulation as `tidal_initial_pos` / `tidal_initial_vel` / `df_traj`.

```julia
using AstroNbodySim, WaveDM
using DataFrames, Unitful, PhysicalParticles
using StructArrays

# 1. Load the catalogue and pick Crater II
df = AstroIC.load_data_MW_satellites()
id = findfirst(df.Galaxy .== "Crater II")

# 2. Convert observables → Galactocentric initial state
pos = PVector(df.X[id], df.Y[id], df.Z[id], u"kpc")
vel = -PVector(df.v_X[id], df.v_Y[id], df.v_Z[id], u"km/s")

# 3. Build a one-particle N-body system
Crater_II = StructArray([Star(uAstro; id = 1)])
Crater_II.Pos[1] = pos
Crater_II.Vel[1] = vel

# 4. Backward integration in a (time-dependent) MW potential
function MW_bg_force(particle, t)
    acc_DM = acc_from_MW_halo(particle)         # user-defined MW DM acceleration
    acc_b  = acc_from_MW_baryons(particle)      # user-defined MW disk+bulge
    return acc_DM + acc_b
end

sim_traj = Simulation(deepcopy(Crater_II);
    bgforce   = Function[MW_bg_force],
    TimeEnd   = 6.0u"Gyr",
    TimeStep  = 0.0u"Gyr",                       # adaptive
)
run(sim_traj)

# 5. The orbit is in `sim_traj.saveresult` — pass it as df_traj to a wave-DM run.
```

The full MW background force should include the disk, bulge, and dark-matter halo (Zhu et al. 2023 form). See `examples/03_lookback.jl` in the repository for a complete working version with the LMC perturbation.

## 5. Self-interacting wave DM

Set $\kappa \neq 0$ to add a contact self-interaction $V = \kappa|\psi|^2$. WaveDM.jl supports repulsive ($\kappa > 0$) and attractive ($\kappa < 0$) regimes.

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model = :MW,
    IC    = (x, y, z) -> exp(-(x^2 + y^2 + z^2) / 0.1^2),
    V     = (x, y, z, ψ) -> 0.0,
    κ     = 50.0,                               # repulsive self-interaction
    Nx = Ny = Nz = 128,
    Xmax = 2.0,
    Tmax = 0.05,
    autoset_timestep = true,
    gpu = true,
    title = "SelfInteracting",
    outputdir = "results/SelfInteracting",
)
```

For $\kappa > 0$ the auto-timestep includes the term $2\pi / (|\kappa|\max|\psi|^2)$ which limits phase variations from the self-interaction.

## 6. Custom potential — cross-disciplinary use

The pluggable `V(x, y, z, ψ)` interface lets WaveDM.jl solve any local nonlinear Schrödinger problem. For example, a 1D effective optical lattice:

```julia
using WaveDM

V_optical(x, y, z, ψ) = 1000 * sin(2π * x)^2 + 1000 * sin(2π * y)^2

simulate_waveDM(;
    model = :MW,
    IC    = (x, y, z) -> exp(-(x^2 + y^2 + z^2) / 0.05^2) * exp(im * 50 * x),
    V     = V_optical,
    κ     = 5000.0,                             # strong Kerr nonlinearity
    Nx = Ny = Nz = 128,
    Xmax = 1.0,
    Tmax = 0.01,
    autoset_timestep = true,
    gpu = true,
    title = "OpticalLattice",
    outputdir = "results/OpticalLattice",
)
```

This makes WaveDM.jl directly usable for nonlinear-optics and BEC problems without modifying the source.

## 7. Multi-level parallel execution

### 7.1 Multi-threaded FFTW on a CPU node

```bash
julia -t 16 my_script.jl
```

```julia
using FFTW
using WaveDM

FFTW.set_provider!("mkl")                      # or "fftw" — must precede the first FFT
FFTW.set_num_threads(Threads.nthreads())

simulate_waveDM(; Nx = Ny = Nz = 384, gpu = false, ...)  # no extra keyword needed
```

### 7.2 Distributed run across multiple cores / nodes

```julia
using Distributed
addprocs(8)
@everywhere using WaveDM

simulate_waveDM(;
    Nx = Ny = Nz = 256,
    Xmax = 5.0,
    Tmax = 2.0u"Gyr",
    autoset_timestep = true,
    distributed_memory = true,                 # use DArray
    gpu = false,
)
```

For multi-node runs, replace `addprocs(8)` with a cluster manager, e.g.

```julia
using Distributed
addprocs([("node1.example.com", :auto), ("node2.example.com", :auto)])
@everywhere using WaveDM
```

### 7.3 GPU run

```julia
using WaveDM
using CUDA
@assert CUDA.functional()

simulate_waveDM(;
    Nx = Ny = Nz = 512,
    Xmax = 5.0,
    Tmax = 5.0u"Gyr",
    autoset_timestep = true,
    gpu = true,
    Realtime = true,
)
```

## 8. Re-using initial conditions across runs

`simulate_waveDM(; IC = path)` accepts a JLD2 path that was written by a previous `save_IC = true` run. This is the recommended way to scan parameters (resolution, $\zeta_{\rm rot}$, etc.) on a fixed initial condition.

```julia
# 1. First run: build and save the IC
simulate_waveDM(;
    model = :MW,
    Xmax = 2.0, Tmax = 0.01, Nt = 1,           # one step is enough if we only want the IC
    Nx = Ny = Nz = 256,
    save_IC = true,
    outputdir = "results/IC_run",
)
# output: results/IC_run/<title>, <suffix> - IC.jld2

# 2. Second run: re-use the saved IC, vary parameters
simulate_waveDM(;
    IC = "results/IC_run/<title>, <suffix> - IC.jld2",
    Xmax = 2.0, Tmax = 0.5,
    autoset_timestep = true,
    rotational_ratio = 0.0,                    # no rotation this time
    bulk_perturb    = false,                   # disable pooling
    gpu = true,
    title = "Scan_rot0",
    outputdir = "results/Scan_rot0",
)
```

## 9. Post-processing — radial profile, virial energies, mass-fraction radii

`simulate_waveDM` returns a `DataFrame dfProp` with the time series of total mass, $R_{10\%}, \ldots, R_{90\%}$ mass-fraction radii, and (if `plot_virial = true`) the kinetic, potential, and quantum energies.

```julia
using CSV, DataFrames, GLMakie

ψ, fig, figMOND, chi2RC, dfProp, dfAcc, avg = simulate_waveDM(;
    Nx = Ny = Nz = 256, plot_virial = true, gpu = true,
)

# 9.1 Mass-fraction radii vs time
fig_r = GLMakie.Figure()
ax = GLMakie.Axis(fig_r[1, 1]; ylabel = "R (kpc)", xlabel = "t (Gyr)")
for i in 1:9
    GLMakie.lines!(ax, dfProp.t, dfProp[!, "R$(i*10)_kpc"]; label = "$(i*10)%")
end
GLMakie.axislegend(ax)
GLMakie.save("radii_vs_t.png", fig_r)

# 9.2 Virial energy convergence
fig_v = GLMakie.Figure()
ax2 = GLMakie.Axis(fig_v[1, 1]; ylabel = "Energy", xlabel = "t (Gyr)")
GLMakie.lines!(ax2, dfProp.t, dfProp.K_total; label = "K")
GLMakie.lines!(ax2, dfProp.t, dfProp.V_total; label = "V")
GLMakie.lines!(ax2, dfProp.t, dfProp.Q_total; label = "Q")
GLMakie.lines!(ax2, dfProp.t, dfProp.E_total; label = "E = K+V+Q")
GLMakie.axislegend(ax2)
GLMakie.save("virial.png", fig_v)
```

A halo is **virialized** when the total energy $\mathcal E = \mathcal K + \mathcal V + \mathcal Q$ is approximately constant and the radii of fixed mass fractions stop growing. The Crater II simulation reaches this state at $t \approx 1.5$ Gyr.

## 10. Best practices

1. **Start small.** Validate with a $64^3$ CPU run before launching a $512^3$ GPU production run.
2. **Enable `autoset_timestep = true`.** It picks $\Delta t$ from the CFL criterion with a 10 % safety margin.
3. **Save initial conditions.** Use `save_IC = true` and re-load with `IC = "..."` to scan parameters on a fixed halo.
4. **Tune `absorb_coeff`.** Defaults of 10 are good for production. Set `absorb_coeff = 0` to disable and see the boundary.
5. **Pick the right `baryon_mode`.** `:mesh` is fastest, `:particles_static` is more flexible, `:particles_dynamic` is the most realistic but most expensive.
6. **Multi-level parallelism is orthogonal.** You can combine `gpu = true` and `distributed_memory = true` for hybrid runs, but in practice most users pick one.
7. **Profile the visualization overhead.** If `Realtime = true` is the bottleneck, switch to `unicode_plot = true` for headless or `Realtime = false` for batch runs.
8. **Check virialization before analyzing.** Use the QSS criteria in [Algorithms §6](algorithms.md#6-quasi-stationary-equilibrium-qss-diagnostics) to verify the halo is in quasi-stationary equilibrium before you fit a profile or compare with observations.
9. **Document your parameters.** Keep a record of `Nx`, `Xmax`, `Tmax`, `mₐ`, `ζ_rot`, `ζ_pool`, `absorb_coeff`, and `baryon_mode` for every run.

## Next steps

- The full list of keyword arguments is in [Simulation API](@ref).
- The mathematical background of the KDK solver is in [Algorithms](@ref).
- If you build something useful on top of WaveDM.jl, consider opening a PR — see the [contribution guide](https://github.com/JuliaAstroSim/WaveDM.jl/blob/main/CONTRIBUTING.md).
