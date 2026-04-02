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
Device configuration for GPU/CPU execution.

# Fields
- `gpu`: Boolean flag indicating whether to use GPU
- `DeviceArray`: Function to convert arrays to device-specific arrays (e.g., cu for GPU)
- `DA`: Function to create distributed arrays (e.g., DArray for distributed memory)
"""
struct DeviceConfig{F1<:Function, F2<:Function}
    gpu::Bool
    DeviceArray::F1
    DA::F2
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
