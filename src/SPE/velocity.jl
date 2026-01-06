# SPE velocity field generation module

# Compute acceleration field from potential
function compute_acceleration_field(ρ_halo, Δ, boundary, gpu, xxx, yyy, zzz, unit_cell_volumn, mass_astro, SofteningLength, potential_astro, pids, length_astro)
    Nx, Ny, Nz = size(ρ_halo)
    if boundary isa Vacuum
        # Use octree to compute gravitational force
        mesh_particles = StructArray(Star(uAstro; id = i+(j-1)*Nx+(k-1)*Nx*Ny) for i in 1:Nx, j in 1:Ny, k in 1:Nz)
        mesh_particles.Pos .= PVector.(xxx * length_astro, yyy * length_astro, zzz * length_astro)
        mesh_particles.Mass .= ρ_halo * unit_cell_volumn * mass_astro

        sim_mesh_force = Simulation(mesh_particles;
            GravitySolver = Tree(),
            pids,
        )

        Φ_WaveDM = compute_potential(sim_mesh_force, mesh_particles.Pos, SofteningLength, Tree(), CPU()) ./ potential_astro
    else
        Φ_WaveDM = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ρ_halo, Periodic(), gpu ? GPU() : CPU()))
    end

    ax_WaveDM, ay_WaveDM, az_WaveDM = grad_central(-Δ..., Φ_WaveDM)

    return Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM
end

# Generate velocity field based on acceleration
function generate_velocity_field(Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM, Φ_b, ax_b, ay_b, az_b, baryon_mode, xxx, yyy, zzz, Δ, velocity_falling, rotational_ratio, velocity_ratio)
    Nx, Ny, Nz = size(xxx)
    if baryon_mode == :ignored
        Φ_all = Φ_WaveDM
        ax_all = ax_WaveDM
        ay_all = ay_WaveDM
        az_all = az_WaveDM
    else
        Φ_all = Φ_WaveDM + Φ_b
        ax_all = ax_WaveDM + ax_b
        ay_all = ay_WaveDM + ay_b
        az_all = az_WaveDM + az_b
    end

    a_CDM = zeros(Nx, Ny, Nz)
    a_CDM[2:end-1,2:end-1,2:end-1] = sqrt.(ax_all[2:end-1,2:end-1,2:end-1].^2 .+ ay_all[2:end-1,2:end-1,2:end-1].^2 .+ az_all[2:end-1,2:end-1,2:end-1].^2)            
    a_CDM[div(end,2),div(end,2),div(end,2)] *= 0

    if velocity_falling
        v = AstroIC.freefall_velocity_acc.(xxx, yyy, zzz, a_CDM) .* velocity_ratio |> collect
    else
        v = AstroIC.rotational_velocity_acc.(xxx, yyy, zzz, a_CDM, rotational_ratio) .* velocity_ratio |> collect
    end

    return v, Φ_all, ax_all, ay_all, az_all
end

# Adjust velocity field (center velocity, bulk perturbation, etc.)
function adjust_velocity_field(v, ρ_halo, bulk_perturb, bulk_size, bulk_shift_size, bulk_center_size)
    Nx, Ny, Nz = size(vx)
    vx = getproperty.(v, :x)
    vy = getproperty.(v, :y)
    vz = getproperty.(v, :z)

    if bulk_perturb
        @info "Setting bulk velocities with size $(bulk_size), shift $(bulk_shift_size)"
        vx = zoom(shrink(vx, bulk_size, bulk_shift_size), bulk_size)
        vy = zoom(shrink(vy, bulk_size, bulk_shift_size), bulk_size)
        vz = zoom(shrink(vz, bulk_size, bulk_shift_size), bulk_size)
    end

    if bulk_center_size > 0
        @info "Setting center velocity to zero, bulk_center_size = $(bulk_center_size)"
        range_x_center = div(Nx,2)-div(bulk_center_size,2)+1:div(Nx,2)-div(bulk_center_size,2)+bulk_center_size
        range_y_center = div(Ny,2)-div(bulk_center_size,2)+1:div(Ny,2)-div(bulk_center_size,2)+bulk_center_size
        range_z_center = div(Nz,2)-div(bulk_center_size,2)+1:div(Nz,2)-div(bulk_center_size,2)+bulk_center_size
        vx[range_x_center, range_y_center, range_z_center] .= 0
        vy[range_x_center, range_y_center, range_z_center] .= 0
        vz[range_x_center, range_y_center, range_z_center] .= 0
    end

    @info "Cancel out the velocity of mass center"
    sum_ρ_halo = sum(ρ_halo)
    vx0 = sum(vx .* ρ_halo) ./ sum_ρ_halo
    vy0 = sum(vy .* ρ_halo) ./ sum_ρ_halo
    vz0 = sum(vz .* ρ_halo) ./ sum_ρ_halo

    vx = vx .- vx0
    vy = vy .- vy0
    vz = vz .- vz0

    return vx, vy, vz
end
