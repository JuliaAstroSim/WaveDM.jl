# Configs

```@meta
CurrentModule = WaveDM
```

WaveDM.jl groups its simulation parameters into a small set of immutable configuration structures. The structures are not strictly required — every field is also exposed as a keyword argument to [`simulate_waveDM`](@ref) / [`SPE3D_waveDM`](@ref) — but they make it easy to thread a complete configuration through helper functions without a long argument list.

All fourteen structures below are exported from the `WaveDM` module.

## 1. Common structures

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    SimulationGrid,
    DeviceConfig,
    TimeStepConfig,
    AstroUnitsConfig,
]
```

### `SimulationGrid`

Carries the grid geometry and the precomputed coordinate arrays that the rest of the code uses. Fields: `Xmax`, `Ymax`, `Zmax`, `Nx`, `Ny`, `Nz`, `Δ`, `x`, `y`, `z`, `xxx`, `yyy`, `zzz`, `r`, `oneMatrix`, `unit_cell_volumn`.

### `DeviceConfig`

Thin compatibility wrapper around [`PhysicalMeshes.ParallelBackend`](https://github.com/JuliaAstroSim/PhysicalMeshes.jl). The preferred way to construct a backend is now [`select_backend`](@ref), which auto-detects the best kind (`:serial` / `:threads` / `:distributed` / `:gpu`) for the current process.

- `backend::ParallelBackend` — the actual backend object.
- `DeviceArray` — convenience alias of `backend.device_array` (the conversion function: `cu` on GPU, `collect` on CPU).
- `DA` — convenience alias of `backend.distribute_array` (the array constructor: `DArray` for distributed runs, `collect` otherwise).
- `gpu::Bool` — `true` dispatches to `CUDA.jl`, `false` stays on CPU.
- `kind::Symbol` — one of `:serial`, `:threads`, `:distributed`, `:gpu`.
- `nthreads`, `nworkers`, `pids` — process topology detected by [`detect_resources`](@ref).

Constructors:

- `DeviceConfig(; gpu=false, distributed=false, kind=:auto, ...)` — full auto-detection, accepting the same kwargs as `select_backend`.
- `DeviceConfig(gpu::Bool, device_array, distribute_array)` — backward-compat positional constructor that builds a minimal backend without auto-detection.

All call sites of the old `.gpu` / `.DeviceArray` / `.DA` fields keep working through property forwarding.

### `TimeStepConfig`

Holds the time vector and step. Fields: `t::Vector{T}`, `dt::T`, `Nt::Int`, `KDK_flag::Bool`.

### `AstroUnitsConfig`

Holds the conversion factors between dimensionless simulation units and physical astrophysical units. Fields: `length_astro`, `time_astro`, `mass_astro`, `density_astro`, `acc_astro`, `velocity_astro`, `potential_astro`, `uT`, `uL`, `uVel`, `uAcc`, `uRho`, `uMomentum`, `h_astro`, `aₛ_astro`, `mₐ_astro`, `c_astro`, `κ_astro`, `G0`, `a0`.

## 2. Kick-step structures

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    GravityConfig,
    TidalFieldConfig,
]
```

### `GravityConfig`

Passed to [`apply_kick_step!`](@ref) each iteration. Fields: `boundary`, `Φ_b`, `sim_mesh_force`, `mesh_particles`, `SofteningLength`, `baryon_mode`, `GravitySolver`, `mass_astro`.

### `TidalFieldConfig`

Bundles the time-dependent MW and LMC tidal pieces. Fields: `MW_tidal_field`, `MW_tidal_interpolate`, `LMC_tidal_field`, `uT`, `tidal_lookback_time`, `df_traj`, `df_traj_LMC`, `length_astro`, `uL`, `MW_grid`, `MW_Phi`, `spl_pot`, `sim_force_baryon`, `particles_LMC`, `sim_traj_LMC`.

## 3. Initial-condition structures

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    InitialConditionsConfig,
    DensityProfileConfig,
    MassRadiusConfig,
]
```

### `InitialConditionsConfig`

Drives [`generate_initial_conditions`](@ref). Fields: `model`, `baryon_mode`, `Np`, `pids`, `bulk_perturb`, `bulk_size`, `bulk_shift_size`, `bulk_center_size`, `reset_velocity`, `static`, `FDM_mass_ratio`, `FDM_radius_ratio`, `rotational_ratio`, `velocity_ratio`, `rotational_ratio_baryon`, `velocity_ratio_baryon`, `velocity_falling`, `MW_disk_RC`, `GravitySolver`, `SofteningLength`.

### `DensityProfileConfig`

The full set of baryonic and halo profile parameters. Fields: `baryon_β`, `baryon_ρ0`, `baryon_r0`, `halo_β`, `halo_ρ0`, `halo_r0`, `halo_α`, `halo_γ`, `halo_Q`, `stellar_TotalMass`, `stellar_ScaleRadius`, `thickness_ratio_stellar`, `gases_TotalMass`, `gases_ScaleRadius`, `thickness_ratio_gases`.

### `MassRadiusConfig`

Simple radial cutoffs used by MOND / RAR diagnostics. Fields: `massRadius`, `minR`, `maxR`.

## 4. Visualization structures

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    VisualizationConfig,
    VisualizationData,
]
```

### `VisualizationConfig`

Controls the live `GLMakie` figure. Fields: `title`, `suffix`, `size`, `StepsBetweenSnapshots`, `Realtime`, `dynamic_colorrange`, `plot_virial`, `plot_optical`, `plot_wavedm`.

### `VisualizationData`

Pre-allocated observables for the per-step update. Fields: `rho`, `rho_max_id`, `total_halo_mass`, `radii`, `r_mass_center`, `target_profile_model`, `target_profile_error`.

## 5. Best-fit extraction structures

```@autodocs
Modules = [WaveDM]
Pages = ["core/configs.jl"]
Order = [:type]
Filter = t -> t in [
    ProfileFitConfig,
    RCFitConfig,
    BestFitConfig,
]
```

### `ProfileFitConfig`

Holds the target density profile (gNFW / Zhao / NFW) and its 1σ uncertainties. Fields: `target_profile_ρ0`, `target_profile_ρ0_u`, `target_profile_ρ0_d`, `target_profile_rs`, `target_profile_rs_u`, `target_profile_rs_d`, `target_profile_α`, `target_profile_α_u`, `target_profile_α_d`, `target_profile_β`, `target_profile_β_u`, `target_profile_β_d`, `target_profile_γ`, `target_profile_γ_u`, `target_profile_γ_d`, `target_fitting_rs_ratio`, `uniform_interval`.

### `RCFitConfig`

Rotation-curve fit and inner-slope $\beta^*$ configuration. Fields: `target_beta_star`, `target_beta_star_u`, `target_beta_star_d`, `target_beta_star_r_min`, `target_beta_star_r_max`, `beta_star_error_threshold`, `Galaxy_i`.

### `BestFitConfig`

Time-averaging and extraction window. Fields: `extract_dwarf_granule`, `extract_min_t`, `extract_mode`, `folder_data`, `average`, `average_start_t`, `average_snapshots`, `average_all`, `average_all_start_t`.

## 6. Usage pattern

The structures are constructed implicitly by `simulate_waveDM` from its keyword arguments. If you are calling `SPE3D_waveDM` directly you can either pass the same keyword arguments or build the structures by hand. For example:

```julia
using WaveDM
using Unitful

config_units = AstroUnitsConfig(
    uconvert(u"kpc",  (8π*(C.h/2π)^2 / (3*mₐ^2 * C.H^2 * Ωₘ₀))^0.25),
    uconvert(u"Gyr",  (3*C.H^2*Ωₘ₀ / (8π))^-0.5),
    # ...
)

config_device = DeviceConfig(true, cu, collect)         # GPU
config_tidal  = TidalFieldConfig(false, true, false,    # no MW / LMC tidal
                                 1.0, 0.0u"Gyr",
                                 DataFrame(), DataFrame(),
                                 1.0, 1.0,
                                 nothing, nothing, nothing,
                                 nothing, nothing, nothing)

# ... then pass these into a custom KDK driver.
```

The keyword-argument interface is the recommended path for users; the structs are most useful for advanced callers who want to encapsulate a "campaign" of related runs in a single object.
