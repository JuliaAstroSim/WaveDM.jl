# Core data structures for WaveDM.jl

# ==============================================================================
# Common structures
# ==============================================================================

"""
Grid and coordinate configuration for simulations.

# Fields
- `Xmax`, `Ymax`, `Zmax`: Maximum grid coordinates in each dimension
- `Nx`, `Ny`, `Nz`: Number of grid points in each dimension
- `Δ`: Grid spacing in each dimension
- `x`, `y`, `z`: 1D coordinate arrays
- `xxx`, `yyy`, `zzz`: 3D coordinate arrays for each dimension
- `r`: 3D radial distance array
- `oneMatrix`: 3D matrix of ones with size (Nx, Ny, Nz)
- `unit_cell_volumn`: Volume of a single grid cell
"""
struct SimulationGrid{T} 
    Xmax::T
    Ymax::T
    Zmax::T
    Nx::Int
    Ny::Int
    Nz::Int
    Δ::Vector{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    xxx::Array{T, 3}
    yyy::Array{T, 3}
    zzz::Array{T, 3}
    r::Array{T, 3}
    oneMatrix::Array{T, 3}
    unit_cell_volumn::T
end

"""
Thin compatibility wrapper around `ParallelBackend` so that the
existing call sites (`config_device.gpu`, `config_device.DeviceArray`,
`config_device.DA`) keep working.
The preferred way to construct a
backend is now [`select_backend`](@ref), which picks the best kind
(`:serial` / `:threads` / `:distributed` / `:gpu`) automatically.

# Fields
- `backend::ParallelBackend`  The actual backend object.
- `DeviceArray::Function`     Convenience alias of `backend.device_array`
- `DA::Function`              Convenience alias of `backend.distribute_array`

Constructors:
- `DeviceConfig(; gpu=false, distributed=false, kind=:auto, ...)`
    - Build via auto-detection.  Accepts the same kwargs as `select_backend`.
- `DeviceConfig(gpu::Bool, device_array, distribute_array)`
    - Backward-compat: build a minimal serial/GPU-style config without
      invoking auto-detection.  The new code should prefer the keyword
      constructor.
"""
struct DeviceConfig
    backend::ParallelBackend
    DeviceArray::Function
    DA::Function
end

# Keyword constructor: full auto-detection
function DeviceConfig(; gpu::Bool=false, distributed::Bool=false, kind::Symbol=:auto, kw...)
    backend = select_backend(; gpu=gpu, distributed=distributed, kind=kind, kw...)
    return DeviceConfig(backend, backend.device_array, backend.distribute_array)
end

# Positional constructor for backward compatibility: build a minimal
# backend of the corresponding kind, no auto-detection.
function DeviceConfig(gpu::Bool, device_array::Function, distribute_array::Function)
    kind = if gpu
        :gpu
    elseif distribute_array === identity
        :serial
    else
        :distributed
    end
    backend = select_backend(; kind=kind, gpu=gpu)
    # Override the auto-detected function fields with the ones the caller
    # passed in.  This keeps old code paths working when the user hands
    # us custom `cu` / `DArray` functions.
    backend = ParallelBackend(
        kind, gpu,
        backend.nthreads, backend.nworkers, backend.pids,
        backend.has_darrays, backend.has_parallelops,
        device_array, distribute_array, backend.local_array, backend.release!,
        backend.fft_threads,
    )
    return DeviceConfig(backend, device_array, distribute_array)
end

# Property forwarding: old call sites read `config_device.gpu` etc.
function Base.getproperty(dc::DeviceConfig, sym::Symbol)
    b = getfield(dc, :backend)
    sym === :gpu              && return b.gpu
    sym === :DeviceArray      && return getfield(dc, :DeviceArray)
    sym === :DA               && return getfield(dc, :DA)
    sym === :device_array     && return b.device_array
    sym === :distribute_array && return b.distribute_array
    sym === :local_array      && return b.local_array
    sym === :release!         && return b.release!
    sym === :kind             && return b.kind
    sym === :nthreads         && return b.nthreads
    sym === :nworkers         && return b.nworkers
    sym === :pids             && return b.pids
    sym === :fft_threads      && return b.fft_threads
    sym === :has_darrays      && return b.has_darrays
    sym === :has_parallelops  && return b.has_parallelops
    sym === :backend          && return b
    return getproperty(b, sym)
end

"""
Time step configuration.
"""
struct TimeStepConfig{T}
    t::Vector{T}
    dt::T
    Nt::Int
    KDK_flag::Bool
end

"""
Astrophysical units configuration.
"""
struct AstroUnitsConfig{T}
    length_astro
    time_astro
    mass_astro
    density_astro
    acc_astro
    velocity_astro
    potential_astro
    uT::T
    uL::T
    uVel::T
    uAcc::T
    uRho::T
    uMomentum::T
    h_astro
    aₛ_astro
    mₐ_astro
    c_astro
    κ_astro
    G0
    a0
end

# ==============================================================================
# Structures for apply_kick_step!
# ==============================================================================

"""
Gravitational configuration for kick step.
"""
# Make GravityConfig more flexible by removing type parameter T
struct GravityConfig
    boundary
    Φ_b
    sim_mesh_force
    mesh_particles
    SofteningLength
    baryon_mode
    GravitySolver
    mass_astro
end

"""
Tidal field configuration for kick step.
"""
struct TidalFieldConfig
    MW_tidal_field::Bool
    MW_tidal_interpolate::Bool
    LMC_tidal_field::Bool
    uT
    tidal_lookback_time
    df_traj
    df_traj_LMC
    length_astro
    uL
    MW_grid
    MW_Phi
    spl_pot
    sim_force_baryon
    particles_LMC
    sim_traj_LMC
end



# ==============================================================================
# Structures for generate_initial_conditions
# ==============================================================================

"""
Initial conditions configuration.
"""
struct InitialConditionsConfig{T}
    model::Symbol
    baryon_mode::Symbol
    Np::Int
    pids
    bulk_perturb::Bool
    bulk_size::Int
    bulk_shift_size::Int
    bulk_center_size::Int
    reset_velocity::Bool
    static::Bool
    FDM_mass_ratio::T
    FDM_radius_ratio::T
    rotational_ratio::T
    velocity_ratio::T
    rotational_ratio_baryon::T
    velocity_ratio_baryon::T
    velocity_falling::Bool
    MW_disk_RC::Bool
    GravitySolver
    SofteningLength
end

"""
Density profile configuration.
"""
# Make DensityProfileConfig more flexible by allowing different types for different fields
struct DensityProfileConfig
    baryon_β
    baryon_ρ0
    baryon_r0
    halo_β
    halo_ρ0
    halo_r0
    halo_α
    halo_γ
    halo_Q
    stellar_TotalMass
    stellar_ScaleRadius
    thickness_ratio_stellar
    gases_TotalMass
    gases_ScaleRadius
    thickness_ratio_gases
end

"""
Mass radius configuration.
"""
struct MassRadiusConfig{T}
    massRadius::T
    minR::T
    maxR::T
end

# ==============================================================================
# Structures for setup_visualization
# ==============================================================================

"""
Visualization configuration.
"""
struct VisualizationConfig
    title::String
    suffix::String
    size::Tuple{Int, Int}
    StepsBetweenSnapshots::Int
    Realtime::Bool
    dynamic_colorrange::Bool
    plot_virial::Bool
    plot_optical::Bool
    plot_wavedm::Bool
end

"""
Visualization data configuration.
"""
struct VisualizationData{T}
    rho::Array{T, 3}
    rho_max_id::CartesianIndex{3}
    total_halo_mass::T
    radii::Vector{T}
    r_mass_center::Array{T, 3}
    target_profile_model::Symbol
    target_profile_error::Bool
end

# ==============================================================================
# Structures for best fit extraction
# ==============================================================================

"""
Profile fitting configuration.
"""
struct ProfileFitConfig{DEN,LEN}
    target_profile_ρ0::DEN
    target_profile_ρ0_u::DEN
    target_profile_ρ0_d::DEN
    target_profile_rs::LEN
    target_profile_rs_u::LEN
    target_profile_rs_d::LEN
    target_profile_α::Real
    target_profile_α_u::Real
    target_profile_α_d::Real
    target_profile_β::Real
    target_profile_β_u::Real
    target_profile_β_d::Real
    target_profile_γ::Real
    target_profile_γ_u::Real
    target_profile_γ_d::Real
    target_fitting_rs_ratio::Real
    uniform_interval::Bool
end

"""
Rotation curve fitting configuration.
"""
struct RCFitConfig{LEN}
    target_beta_star::Real
    target_beta_star_u::Real
    target_beta_star_d::Real
    target_beta_star_r_min::LEN
    target_beta_star_r_max::LEN
    beta_star_error_threshold::Real
    Galaxy_i::Int
end

"""
Best fit extraction configuration.
"""
struct BestFitConfig{T}
    extract_dwarf_granule::Bool
    extract_min_t::T
    extract_mode::Symbol
    folder_data::String
    average::Bool
    average_start_t::T
    average_snapshots::Bool
    average_all::Bool
    average_all_start_t::T
end
