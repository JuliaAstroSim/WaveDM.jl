"""
$(TYPEDSIGNATURES)

Apply a full kick step including gravitational potential computation.
This function encapsulates the entire kick step from density computation to potential application.

# Arguments
- `device_ψ`: Device-resident wavefunction.  *Consumed* by this function:
  the underlying GPU buffer (if any) is freed before the function returns and
  the Julia binding on the caller side is invalidated.
- `ψ`: Host-resident wavefunction (used for evaluating the external potential `V`).
- `V`: Potential function with signature `V(x, y, z, ψ) -> scalar`.
- `rho_max`: Maximum density value (kept for backward compatibility; not used).
- `rho_max_id`: Cartesian index of maximum density (used as the centering point for the tidal potential correction).
- `grid`: SimulationGrid (contains: xxx, yyy, zzz, Δ, Nx, Ny, Nz, unit_cell_volumn, oneMatrix).
- `config_gravity`: GravityConfig (constants only).
- `config_tidal`: TidalFieldConfig (constants only).
- `config_device`: DeviceConfig (constants only: gpu, DeviceArray).
- `dt`: Time step (variable, passed as parameter).
- `i`: Current iteration index (variable, passed as parameter).
- `t`: Time vector (variable, passed as parameter).

# Returns
- `Φ_all`: Host-resident total gravitational potential.
- `Φ_WaveDM`: Host-resident WaveDM-only gravitational potential.
- `spec`: Device-resident Fourier transform of the kicked wavefunction, ready for the drift step.

# Memory management
On the GPU path, the temporary device buffers that this function creates
(`device_Φ_all`, `device_Φ_WaveDM`, `potential_grav_static`, `device_V`,
`potential_all`, `nonlinear_term`, and `device_ψ`) are released through
`CUDA.unsafe_free!` as soon as they are no longer needed so that peak
GPU memory stays close to the size of one wavefunction array.  On the CPU
path Julia's GC will reclaim the corresponding host arrays when their
bindings go out of scope.
"""
function apply_kick_step!(
    device_ψ,
    ψ,
    V::F,
    rho,
    rho_max,
    rho_max_id::CartesianIndex{3},
    grid::SimulationGrid,
    config_gravity::GravityConfig,
    config_tidal::TidalFieldConfig,
    config_device::DeviceConfig,
    dt::Real,
    i::Int,
    t::Vector{<:Real}
) where {F<:Function}
    device_Φ_all, device_Φ_WaveDM = compute_gravitational_potential(device_ψ, config_gravity, grid, config_device)

    Φ_all = collect(device_Φ_all)
    Φ_WaveDM = collect(device_Φ_WaveDM)
    _release!(device_Φ_all, config_device.gpu)
    _release!(device_Φ_WaveDM, config_device.gpu)

    if config_tidal.MW_tidal_field
        add_tidal_potential!(Φ_all, rho, config_tidal, grid, config_gravity, i, t, rho_max_id)
    end

    potential_grav_static = config_device.DeviceArray(Φ_all)
    device_V = config_device.DeviceArray(V.(grid.xxx, grid.yyy, grid.zzz, ψ))
    potential_all = device_V .+ potential_grav_static
    _release!(potential_grav_static, config_device.gpu)
    _release!(device_V, config_device.gpu)

    nonlinear_term = exp.(-im .* dt .* potential_all) .* device_ψ
    _release!(device_ψ, config_device.gpu)
    _release!(potential_all, config_device.gpu)

    spec = fft(nonlinear_term)
    _release!(nonlinear_term, config_device.gpu)

    return Φ_all, Φ_WaveDM, spec
end

"""
$(TYPEDSIGNATURES)

Apply a full drift step including boundary condition.
This function encapsulates the entire drift step from Fourier space transformation to boundary application.

`spec` and `linear_phase` are consumed by this function; their bindings on
the caller side should not be reused after the call.
"""
function apply_drift_step!(spec, linear_phase, border, gpu::Bool, DeviceArray)
    device_linear_phase = DeviceArray(linear_phase)
    spec .*= device_linear_phase
    newψ = ifft(spec)
    _release!(spec, gpu)
    _release!(device_linear_phase, gpu)

    device_border = DeviceArray(border)
    device_ψ = device_border .* newψ
    _release!(newψ, gpu)
    _release!(device_border, gpu)

    return device_ψ
end

# Internal helper: release a device buffer when running on the GPU, otherwise
# let Julia's GC handle the host-side array.  Wrapping the conditional in a
# single function avoids the misleading `x = nothing` branches that were
# scattered through the kick / drift steps and clearly documents the intent
# at every call site.
@inline function _release!(buf, gpu::Bool)
    if gpu
        CUDA.unsafe_free!(buf)
    end
    return nothing
end
