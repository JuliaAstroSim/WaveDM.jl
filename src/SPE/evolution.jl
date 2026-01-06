# SPE time evolution module

"""
$(TYPEDSIGNATURES)

Compute gravitational potential from density.
"""
function compute_gravitational_potential(device_ψ, boundary, Δ, Nx, Ny, Nz, unit_cell_volumn, gpu, Φ_b, DeviceArray,
                                         sim_mesh_force, mesh_particles, SofteningLength, potential_astro, baryon_mode, mass_astro, rho_max, rho_max_id)
    ψ² = abs.(device_ψ).^2
    ψ² .-= mean(ψ²)

    if boundary isa Vacuum
        sim_mesh_force.simdata.tree.data.Mass .= ψ²[:] * unit_cell_volumn * mass_astro
        AstroNbodySim.rebuild_tree(sim_mesh_force)
        Φ_WaveDM = compute_potential(sim_mesh_force, mesh_particles.Pos, SofteningLength, Tree(), CPU()) ./ potential_astro
        potential_grav = DeviceArray(Φ_WaveDM) + DeviceArray(Φ_b)
    else
        if baryon_mode == :ignored
            potential_grav = 4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ψ², Periodic(), gpu ? GPU() : CPU()) #TODO: optimize performance
        else
            device_Φ_b = DeviceArray(Φ_b)
            Φ_WaveDM = 4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ψ², Periodic(), gpu ? GPU() : CPU())
            potential_grav = Φ_WaveDM + device_Φ_b
            
            gpu && CUDA.unsafe_free!(device_Φ_b)
            gpu && CUDA.unsafe_free!(Φ_WaveDM)
        end
    end
    
    return potential_grav
end

"""
$(TYPEDSIGNATURES)

Apply a full kick step including gravitational potential computation.
This function encapsulates the entire kick step from density computation to potential application.
"""
function apply_kick_step!(device_ψ, ψ, V, xxx, yyy, zzz, boundary, Δ, Nx, Ny, Nz, MW_grid, MW_Phi, spl_pot, sim_force_baryon,
                              unit_cell_volumn, gpu, Φ_b, DeviceArray, sim_mesh_force, MW_tidal_interpolate, LMC_tidal_field, t, i, uT, tidal_lookback_time, df_traj, df_traj_LMC, length_astro, uL,
                              mesh_particles, SofteningLength, potential_astro, baryon_mode, GravitySolver, particles_LMC, sim_traj_LMC, rho_max_id, oneMatrix,
                              mass_astro, dt)
        potential_grav = compute_gravitational_potential(device_ψ, boundary, Δ, Nx, Ny, Nz, unit_cell_volumn, gpu, Φ_b, DeviceArray, sim_mesh_force, mesh_particles, SofteningLength, potential_astro, baryon_mode, mass_astro, rho_max, rho_max_id)
        
        Φ_all = collect(potential_grav)
        gpu ? CUDA.unsafe_free!(potential_grav) : (potential_grav=nothing)

        MW_tidal_field && add_tidal_potential!(Φ_all, MW_tidal_interpolate, LMC_tidal_field, t, i, uT, tidal_lookback_time, df_traj, df_traj_LMC, xxx, yyy, zzz, length_astro, uL, MW_grid, MW_Phi, spl_pot, sim_force_baryon, SofteningLength, GravitySolver, particles_LMC, sim_traj_LMC, rho_max_id, oneMatrix, Nx, Ny, Nz)            
        flag_cancel_shift && cancel_field_gradient_at_center!(Φ_all, rho_max_id, oneMatrix, Nx, Ny, Nz)
        
        potential_grav_static = DeviceArray(Φ_all)
        device_V = DeviceArray(V.(xxx, yyy, zzz, ψ))
        potential_all = device_V + potential_grav_static
        gpu ? CUDA.unsafe_free!(potential_grav_static) : (potential_grav_static=nothing)
        gpu ? CUDA.unsafe_free!(device_V) : (device_V=nothing) 

        nonlinear_term = exp.(-im*dt*potential_all) .* device_ψ
        gpu ? CUDA.unsafe_free!(device_ψ) : (device_ψ=nothing) 
        gpu ? CUDA.unsafe_free!(potential_all) : (potential_all=nothing) 

        spec = fft(nonlinear_term)
        gpu ? CUDA.unsafe_free!(nonlinear_term) : (nonlinear_term=nothing)
    
    return device_ψ, Φ_all, spec, rho_max, rho_max_id
end

"""
$(TYPEDSIGNATURES)

Apply a full drift step including boundary condition.
This function encapsulates the entire drift step from Fourier space transformation to boundary application.
"""
function apply_drift_step!(spec, linear_phase, boarder, gpu, DeviceArray)
    device_linear_phase = DeviceArray(linear_phase)
    spec .*= device_linear_phase
    newψ = ifft(spec)
    gpu ? CUDA.unsafe_free!(spec) : (spec = nothing)
    gpu ? CUDA.unsafe_free!(device_linear_phase) : (device_linear_phase = nothing)

    device_boarder = DeviceArray(boarder)
    device_ψ = device_boarder .* newψ
    gpu ? CUDA.unsafe_free!(newψ) : (newψ = nothing)
    gpu ? CUDA.unsafe_free!(device_boarder) : (device_boarder = nothing)
    
    return device_ψ
end

# Export functions
export apply_kick_step!, apply_drift_step!