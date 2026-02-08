"""
$(TYPEDSIGNATURES)

Compute gravitational potential from density using config objects.
Accesses config fields directly without internal unpacking.
"""
function compute_gravitational_potential(
    device_ψ,
    config_gravity::GravityConfig,
    grid::SimulationGrid,
    config_device::DeviceConfig,
)
    ψ² = abs.(device_ψ).^2
    ψ² .-= mean(ψ²)

    if config_gravity.boundary isa Vacuum
        config_gravity.sim_mesh_force.simdata.tree.data.Mass .= ψ²[:] * grid.unit_cell_volumn * config_gravity.mass_astro
        AstroNbodySim.rebuild_tree(config_gravity.sim_mesh_force)
        device_Φ_WaveDM = compute_potential(config_gravity.sim_mesh_force, config_gravity.mesh_particles.Pos, config_gravity.SofteningLength, Tree(), CPU()) ./ config_gravity.mass_astro
        device_Φ_all = config_device.DeviceArray(device_Φ_WaveDM) + config_device.DeviceArray(config_gravity.Φ_b)
    else
        if config_gravity.baryon_mode == :ignored
            device_Φ_all = device_Φ_WaveDM = 4π * fft_poisson(grid.Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ψ², Periodic(), config_device.gpu ? GPU() : CPU())
        else
            device_Φ_b = config_device.DeviceArray(config_gravity.Φ_b)
            device_Φ_WaveDM = 4π * fft_poisson(grid.Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ψ², Periodic(), config_device.gpu ? GPU() : CPU())
            device_Φ_all = device_Φ_WaveDM + device_Φ_b
            
            config_device.gpu && CUDA.unsafe_free!(device_Φ_b)
            config_device.gpu && CUDA.unsafe_free!(device_Φ_WaveDM)
        end
    end
    
    return device_Φ_all, device_Φ_WaveDM
end

function baryon_add_WaveDM_acc(sim_force_baryon, config_mesh, meshpos, meshacc_WaveDM)
    baryon_particles_temp = get_local_data(sim_force_baryon) # access the Array pointer
    for k in eachindex(baryon_particles_temp)
        pos = baryon_particles_temp[k].Pos
        acc = uconvert(u"kpc/Gyr^2", baryon_particles_temp[k].Acc + mesh2particle(meshpos, config_mesh, meshacc_WaveDM, pos, config_mesh.mode, config_mesh.assignment))
        setproperty!!(baryon_particles_temp[k], :Acc, acc)
    end
end