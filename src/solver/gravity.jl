"""
$(TYPEDSIGNATURES)

Compute gravitational potential from density using config objects.
Accesses config fields directly without internal unpacking.
"""
function compute_gravitational_potential(
    device_ψ,
    gravity_config::GravityConfig,
    grid::SimulationGrid,
    device_config::DeviceConfig,
)
    ψ² = abs.(device_ψ).^2
    ψ² .-= mean(ψ²)

    if gravity_config.boundary isa Vacuum
        gravity_config.sim_mesh_force.simdata.tree.data.Mass .= ψ²[:] * grid.unit_cell_volumn * gravity_config.mass_astro
        AstroNbodySim.rebuild_tree(gravity_config.sim_mesh_force)
        Φ_WaveDM = compute_potential(gravity_config.sim_mesh_force, gravity_config.mesh_particles.Pos, gravity_config.SofteningLength, Tree(), CPU()) ./ gravity_config.mass_astro
        potential_grav = device_config.DeviceArray(Φ_WaveDM) + device_config.DeviceArray(gravity_config.Φ_b)
    else
        if gravity_config.baryon_mode == :ignored
            potential_grav = 4π * fft_poisson(grid.Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ψ², Periodic(), device_config.gpu ? GPU() : CPU())
        else
            device_Φ_b = device_config.DeviceArray(gravity_config.Φ_b)
            Φ_WaveDM = 4π * fft_poisson(grid.Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ψ², Periodic(), device_config.gpu ? GPU() : CPU())
            potential_grav = Φ_WaveDM + device_Φ_b
            
            device_config.gpu && CUDA.unsafe_free!(device_Φ_b)
            device_config.gpu && CUDA.unsafe_free!(Φ_WaveDM)
        end
    end
    
    return potential_grav
end
