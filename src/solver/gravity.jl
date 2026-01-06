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
