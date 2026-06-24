# Simulation

```@meta
CurrentModule = WaveDM
```

## Main Simulation Functions

### `simulate_waveDM`

```@docs
WaveDM.simulate_waveDM
```

The high-level entry point. It builds the initial conditions (using `AstroIC.jl`), creates a `SimulationGrid`, calls `SPE3D_waveDM` to evolve the SPE, then handles visualization and post-processing. All of the keyword arguments listed below can be supplied directly to `simulate_waveDM`.

### `SPE3D_waveDM`

```@docs
WaveDM.SPE3D_waveDM
```

The lower-level entry point. It expects you to provide the initial wavefunction `IC` and (optionally) the baryonic density/potential/forces yourself. Use `simulate_waveDM` unless you are integrating WaveDM.jl into a larger pipeline.

## Simulation Workflow

```mermaid
flowchart LR
    A[Grid setup] --> B[Initial conditions]
    B --> C[Baryon setup<br/>mesh / static / dynamic]
    C --> D[IC potential]
    D --> E[Visualization init]
    E --> F[KDK main loop]
    F -->|every StepsBetweenSnapshots| E
    F --> G[Outputs<br/>CSV, JLD2, PNG, MP4]
```

1. **Grid setup** — uniform Cartesian grid (`Nx × Ny × Nz`, side lengths `2*Xmax` etc.).
2. **Initial conditions** — `IC` (initial wavefunction) and `V` (external potential) are evaluated on the grid.
3. **Baryon setup** — depending on `baryon_mode`, baryons are placed on the mesh, sampled as static N-body particles, or co-evolved as dynamic N-body particles.
4. **IC potential** — the gravitational potential $\Phi_{\rm total}$ at $t = 0$ is computed (FFT for periodic, James 1977 for vacuum).
5. **Visualization init** — observables for `Realtime` and the per-snapshot update are wired up.
6. **KDK main loop** — `apply_kick_step!` → `apply_drift_step!` for `Nt` steps.
7. **Outputs** — per-step and final files (CSV for properties, JLD2 for $\psi$/$\Phi$, PNG/MP4 for visualization).

## Key Parameters

### Grid and resolution

| Keyword | Default | Description |
| --- | --- | --- |
| `Xmax`, `Ymax`, `Zmax` | `5`, `Xmax`, `Xmax` | Half-side length of the box in dimensionless units (full side `L = 2*Xmax`) |
| `Nx`, `Ny`, `Nz` | `64`, `Nx`, `Nx` | Grid points per axis; powers of two are fastest |
| `Tmax` | `1` | Total simulation time in dimensionless units (or `u"Gyr"` if `Unitful` is loaded) |
| `Nt` | `32` | Number of time steps (overridden when `autoset_timestep = true`) |

### Initial conditions

| Keyword | Default | Description |
| --- | --- | --- |
| `IC` | Gaussian $(x,y,z)\mapsto e^{-x^2-y^2-z^2}$ | Initial wavefunction — a function of $(x,y,z)$ returning a complex number |
| `V` | `(x,y,z,ψ)->0` | External / self-interaction potential $V(\mathbf r, t, \psi)$ |
| `κ` | `0` | Self-interaction strength. `0` (free), `> 0` (repulsive), `< 0` (attractive) |
| `baryon` | `0` | Static baryon density (used only by `baryon_mode = :mesh`) |
| `mₐ` | `0.2×10⁻²² eV/c²` | Boson mass in physical units (drives the dimensionless scale) |
| `Ωₘ₀` | `0.31` | Present-day matter density parameter |
| `MW_disk_RC` | `false` | For `model = :MW`: drive the velocity field from the MW rotation curve |
| `baryon_fraction_limit` | `1` | Upper bound on the baryon-to-total mass ratio used when generating `model = :dwarf` ICs |

### Baryons

| Keyword | Default | Description |
| --- | --- | --- |
| `baryon_mode` | `:mesh` | One of `:ignored`, `:mesh`, `:particles_static`, `:particles_dynamic` |
| `baryon_potential` | `nothing` | Pre-computed $\Phi_{\rm baryon}$; skip the Poisson solve |
| `baryon_particles` | `nothing` | `StructArray` of N-body particles (`Star` / `Gas`) for static or dynamic modes |
| `SofteningLength` | `1.0u"kpc"` | Gravitational softening length for the N-body force |
| `GravitySolver` | `Tree()` | `Tree()` ($\mathcal O(N\log N)$, Gadget-2) or `DirectSum()` ($\mathcal O(N^2)$) |

`baryon_mode` values:

- `:ignored` — no baryons. Default if you are simulating a dark-matter-only halo.
- `:mesh` — baryons live on the same FFT grid; potential computed by the same FFT Poisson solver as the DM. Best for spatially extended baryons (e.g. dSph, ETG).
- `:particles_static` — baryons are sampled once as N-body particles and held fixed; their gravitational potential is added to $\Phi_{\rm total}$ at every step.
- `:particles_dynamic` — baryons are co-evolved with the DM; their potential is re-computed at every step and the DM acceleration is interpolated back to the particles. Most expensive, most realistic.

### Tidal field

| Keyword | Default | Description |
| --- | --- | --- |
| `MW_tidal_field` | `false` | Add a time-dependent MW tidal field to $\Phi_{\rm total}$ |
| `LMC_tidal_field` | `false` | Add LMC perturbation (requires `MW_tidal_field = true`) |
| `MW_tidal_interpolate` | `true` | Use a pre-computed splined MW potential grid (`MW_pot`) for speed |
| `MW_pot` | `nothing` | Path to a JLD2 file or a pre-built potential grid |
| `MW_pot_Xmax`, `MW_pot_Ymax`, `MW_pot_Zmax` | `100.0u"kpc"` | Half-side lengths of the MW potential interpolation box |
| `MW_pot_N` | `512` | Resolution of the MW potential interpolation box |
| `export_MW_pot` | `false` | If `true`, save the MW tidal-field interpolation grid to a JLD2 file |
| `tidal_BG_sim` | `nothing` | Pre-built MW background N-body `Simulation` (bypasses `MW_pot`) |
| `tidal_lookback_time` | `0.0u"Gyr"` | How far back in time the tidal field is applied |
| `tidal_initial_pos`, `tidal_initial_vel` | Crater II defaults | Initial Galactocentric position / velocity of the satellite |
| `particles_LMC` | `nothing` | `StructArray` of LMC particles (alternative to analytic Plummer) |

### Boundary conditions

| Keyword | Default | Description |
| --- | --- | --- |
| `boundary` | `Periodic()` | `Periodic()` (FFT-based) or `Vacuum()` (James 1977, $\mathcal O(N^4)$) |
| `absorb_coeff` | `0` | Absorbing-layer coefficient $\alpha$; `0` disables the absorber |

### Timestep

| Keyword | Default | Description |
| --- | --- | --- |
| `autoset_timestep` | `false` | Choose $\Delta t$ from the CFL criterion automatically |
| `autoset_timestep_ratio` | `0.9` | Safety factor $\zeta_t$ for the auto-timestep |

### Performance

| Keyword | Default | Description |
| --- | --- | --- |
| `gpu` | `true` | Use `CUDA.jl` for FFT and element-wise work |
| `distributed_memory` | `false` | Use `Distributed.jl` `DArray` for grid arrays |
| `pids` | `workers()` | Worker process IDs for distributed runs |
| `KDK_flag` | `true` | If `true`, use the second-order kick–drift–kick scheme |

### Visualization

| Keyword | Default | Description |
| --- | --- | --- |
| `Realtime` | `true` | Open an interactive `GLMakie` window |
| `unicode_plot` | `false` | Render a terminal heat-map of the density slice via `UnicodePlots.jl` |
| `unicode_heatmap_width` | `150` | Width in characters of the terminal heat-map |
| `StepsBetweenSnapshots` | `5` | Update the visualization and save frames every N steps |
| `size` | `(2400, 1400)` | Figure size of the overview plot |
| `plotOptical` | `false` | Plot baryonic (optical) properties in addition to wave-DM |
| `plotWaveDM` | `false` | Plot the wave-DM-only view alongside the total |
| `plot_virial` | `false` | Add a virial-energy axis to the live figure |
| `dynamic_colorrange` | `true` | Auto-scale the colormap to the running maximum density |
| `extract_dwarf_granule` | `false` | Post-process: extract the snapshot best fitting a target profile/RC |
| `extract_mode` | `:profile` | `:profile` or `:RC` extraction mode |
| `extract_min_t` | `0.2u"Gyr"` | Earliest simulation time considered for extraction |
| `best_fit_halo_mass` | `false` | If `true`, optimize the halo mass during the RC fit |

### Output

| Keyword | Default | Description |
| --- | --- | --- |
| `title` | `"WaveDM"` | Run title (used in output filenames) |
| `outputdir` | `joinpath(@__DIR__, "output")` | Where CSV / JLD2 / PNG / MP4 files are written |
| `save_IC` | `true` | Save the initial conditions to disk |
| `save_phi` | `true` | Save the final $\Phi$ and $\psi$ to JLD2 |
| `export_particles` | `false` | If `true`, return the sampled baryon particles (no SPE run) |
| `IC_only` | `false` | If `true`, return the initial wavefunction (no SPE run) |

## Minimal Working Examples

### 1. Smoke test (CPU)

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
    unicode_plot = true,
    outputdir = mktempdir(),
)
```

### 2. GPU production run

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model = :MW,
    Xmax = 5.0,
    Tmax = 5.0u"Gyr",
    Nx = Ny = Nz = 256,
    Nt = 5000,
    autoset_timestep = true,
    gpu = true,
    Realtime = true,
    plot_virial = true,
    title = "MW_256_T5Gyr",
    outputdir = "results/MW_256_T5Gyr",
)
```

### 3. Crater II ultra-faint dwarf

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model       = :dwarf_UFDs,
    Galaxy_id   = 6,                 # Crater II
    V           = (x, y, z, ψ) -> 0.0,
    Nx          = 384,
    Xmax        = 20u"kpc",
    Tmax        = 6.0u"Gyr",
    autoset_timestep = true,
    absorb_coeff = 10,
    gpu         = true,
    Realtime    = true,
    title       = "CraterII",
)
```

### 4. Crater II with MW tidal field

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
    absorb_coeff = 10,
    MW_tidal_field = true,
    LMC_tidal_field = true,         # include LMC perturbation
    gpu         = true,
    title       = "CraterII_MW_LMC",
)
```

### 5. Baryons as dynamic N-body particles

```julia
using WaveDM
using PhysicalParticles, Unitful

baryons = generate(AstroIC.ExponentialDisk(; ...))  # see AstroIC.jl

simulate_waveDM(;
    Nx = 256,
    Xmax = 5.0,
    Tmax = 2.0u"Gyr",
    autoset_timestep = true,
    baryon_mode       = :particles_dynamic,
    baryon_particles  = baryons,
    GravitySolver     = Tree(),
    SofteningLength   = 1.0u"kpc",
    gpu = true,
)
```

### 6. Distributed run

```julia
using Distributed
addprocs(8)
@everywhere using WaveDM

simulate_waveDM(;
    Nx = Ny = Nz = 256,
    Tmax = 2.0u"Gyr",
    autoset_timestep = true,
    distributed_memory = true,         # use DArray for grid arrays
    gpu = false,                       # GPU and distributed are not orthogonal but mixing
                                      # them is an advanced use case
)
```

## Return Values

```julia
ψ, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2 = simulate_waveDM(...)
```

| Return | Description |
| --- | --- |
| `ψ` | Final wavefunction (complex array, dimensionless units) |
| `fig` | Overview Makie figure (snapshot of the live plot at end of run) |
| `figMOND` | MOND / RAR diagnostic figure, or `nothing` if `baryon_mode = :ignored` |
| `chi2RC` | $\chi^2$ of the rotation-curve fit, or `nothing` for DM-only runs |
| `dfProp` | DataFrame with the time series of total mass, mass-fraction radii, virial energies, momenta |
| `dfAcc` | DataFrame with the radial acceleration profile (DM-only or DM+baryon+MOND) |
| `averaged_ψ2` | Time-averaged $|\psi|^2$ if `average = true` was passed, else `nothing` |

## Output Files

The simulation writes several file types into `outputdir`:

- `$(title), $(suffix) - Overview IC.png` — initial snapshot
- `$(title), $(suffix) - Overview.mp4` — animation of the entire run (if `Realtime = true`)
- `$(title), $(suffix) - Overview Prop.png` — final-state overview
- `$(title), $(suffix) - RC.png` / `RAR.png` — rotation curve / radial acceleration relation
- `$(title), $(suffix) - Prop best fit.jld2` — best-fit snapshot (if `extract_dwarf_granule = true`)
- `$(title), $(suffix) - acc.csv` / `acc best fit.csv` — radial acceleration data
- `$(title), $(suffix) - best_fit_time.csv` — best-fit simulation time
- `timestep_*.jld2` — per-snapshot ψ / Φ (if `average_snapshots = true`)

The `suffix` is a human-readable summary of the run, e.g. `Nx(384), Xmax(20.00), Nt(377), Tmax(6.00), DM_m(1.00)`.

## Performance Considerations

- **GPU vs CPU.** `gpu = true` is most useful for grids $\gtrsim 256^3$. For tiny problems the kernel-launch overhead dominates and CPU is faster.
- **Grid size.** Powers of two (64, 128, 256, …) maximize FFT performance. The code also runs on non-power-of-two grids.
- **Timestep.** `autoset_timestep = true` keeps $\Delta t$ near the CFL limit without manual tuning.
- **`StepsBetweenSnapshots`.** Larger values skip more visualization updates per step, which is useful when the plot is the bottleneck.
- **Memory.** On GPU, $N^3$ complex arrays dominate. A $512^3$ single-precision grid is $\sim 2$ GiB; a $1024^3$ grid is $\sim 16$ GiB.
- **`KDK_flag`.** Leave at `true`. Setting it to `false` reduces the scheme to first-order Lie splitting and is provided only for code-experimentation purposes.

## Related Modules

- [Configs](configs.md) — configuration structures
- [KDK](KDK.md) — the kick–drift–kick implementation
