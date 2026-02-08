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
- `config_gravity`: GravityConfig (constants only)
- `config_tidal`: TidalFieldConfig (constants only)
- `config_device`: DeviceConfig (constants only: gpu, DeviceArray)
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
    config_gravity::GravityConfig,
    config_tidal::TidalFieldConfig,
    config_device::DeviceConfig,
    dt::Real,
    i::Int,
    t::Vector{<:Real}
)
    device_Φ_all, device_Φ_WaveDM = compute_gravitational_potential(device_ψ, config_gravity, grid, config_device)
    
    Φ_all = collect(device_Φ_all)
    Φ_WaveDM = collect(device_Φ_WaveDM)
    config_device.gpu ? CUDA.unsafe_free!(device_Φ_all) : (device_Φ_all = nothing)
    config_device.gpu ? CUDA.unsafe_free!(device_Φ_WaveDM) : (device_Φ_WaveDM = nothing)
    
    if config_tidal.MW_tidal_field
        add_tidal_potential!(Φ_all, config_tidal, grid, config_gravity, i, t)
        cancel_field_gradient_at_center!(Φ_all, rho_max_id, grid.oneMatrix, grid.Nx, grid.Ny, grid.Nz)
    end
    
    potential_grav_static = config_device.DeviceArray(Φ_all)
    device_V = config_device.DeviceArray(V.(grid.xxx, grid.yyy, grid.zzz, ψ))
    potential_all = device_V + potential_grav_static
    config_device.gpu ? CUDA.unsafe_free!(potential_grav_static) : (potential_grav_static = nothing)
    config_device.gpu ? CUDA.unsafe_free!(device_V) : (device_V = nothing)
    
    nonlinear_term = exp.(-im * dt * potential_all) .* device_ψ
    config_device.gpu ? CUDA.unsafe_free!(device_ψ) : (device_ψ = nothing)
    config_device.gpu ? CUDA.unsafe_free!(potential_all) : (potential_all = nothing)
    
    spec = fft(nonlinear_term)
    config_device.gpu ? CUDA.unsafe_free!(nonlinear_term) : (nonlinear_term = nothing)
    
    return Φ_all, Φ_WaveDM, spec
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
