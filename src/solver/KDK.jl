"""
$(TYPEDSIGNATURES)

Apply a full kick step including gravitational potential computation.
This function encapsulates the entire kick step from density computation to potential application.

# Arguments
- `device_ψ`: Device-resident wavefunction
- `ψ`: Host-resident wavefunction
- `V`: Potential function
- `rho_max`: Maximum density value
- `rho_max_id`: Cartesian index of maximum density
- `grid`: SimulationGrid (contains: xxx, yyy, zzz, Δ, Nx, Ny, Nz, unit_cell_volumn, oneMatrix)
- `gravity_config`: GravityConfig (constants only)
- `tidal_config`: TidalFieldConfig (constants only)
- `device_config`: DeviceConfig (constants only: gpu, DeviceArray)
- `dt`: Time step (variable, passed as parameter)
- `i`: Current iteration index (variable, passed as parameter)
- `t`: Time vector (variable, passed as parameter)
"""
function apply_kick_step!(
    device_ψ,
    ψ,
    V::Function,
    rho_max,
    rho_max_id::CartesianIndex{3},
    grid::SimulationGrid,
    gravity_config::GravityConfig,
    tidal_config::TidalFieldConfig,
    device_config::DeviceConfig,
    dt::Real,
    i::Int,
    t::Vector{<:Real}
)::Tuple{Array{Float64, 3}, Any}
    potential_grav = compute_gravitational_potential(device_ψ, gravity_config, grid, device_config)
    
    Φ_all = collect(potential_grav)
    device_config.gpu ? CUDA.unsafe_free!(potential_grav) : (potential_grav = nothing)
    
    if tidal_config.MW_tidal_field
        add_tidal_potential!(Φ_all, tidal_config, grid, gravity_config, i, t)
        cancel_field_gradient_at_center!(Φ_all, rho_max_id, grid.oneMatrix, grid.Nx, grid.Ny, grid.Nz)
    end
    
    potential_grav_static = device_config.DeviceArray(Φ_all)
    device_V = device_config.DeviceArray(V.(grid.xxx, grid.yyy, grid.zzz, ψ))
    potential_all = device_V + potential_grav_static
    device_config.gpu ? CUDA.unsafe_free!(potential_grav_static) : (potential_grav_static = nothing)
    device_config.gpu ? CUDA.unsafe_free!(device_V) : (device_V = nothing)
    
    nonlinear_term = exp.(-im * dt * potential_all) .* device_ψ
    device_config.gpu ? CUDA.unsafe_free!(device_ψ) : (device_ψ = nothing)
    device_config.gpu ? CUDA.unsafe_free!(potential_all) : (potential_all = nothing)
    
    spec = fft(nonlinear_term)
    device_config.gpu ? CUDA.unsafe_free!(nonlinear_term) : (nonlinear_term = nothing)
    
    return Φ_all, spec
end

"""
$(TYPEDSIGNATURES)

Apply a full drift step including boundary condition.
This function encapsulates the entire drift step from Fourier space transformation to boundary application.
"""
function apply_drift_step!(spec, linear_phase, border, gpu, DeviceArray)
    device_linear_phase = DeviceArray(linear_phase)
    spec .*= device_linear_phase
    newψ = ifft(spec)
    gpu ? CUDA.unsafe_free!(spec) : (spec = nothing)
    gpu ? CUDA.unsafe_free!(device_linear_phase) : (device_linear_phase = nothing)
    
    device_border = DeviceArray(border)
    device_ψ = device_border .* newψ
    gpu ? CUDA.unsafe_free!(newψ) : (newψ = nothing)
    gpu ? CUDA.unsafe_free!(device_border) : (device_border = nothing)
    
    return device_ψ
end
