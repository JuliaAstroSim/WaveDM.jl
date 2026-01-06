"""
$(TYPEDSIGNATURES)

Compute total Milky Way density at position (x, y, z).
"""
function density_all_MW(x, y, z)
    ρ_bulge, ρ_discs, ρ_gas, ρ_halo = GalacticDynamics.milkyway_zhu2023how()

    function ρ_all(x, y, z)
        R = sqrt(x^2 + y^2)
        r = sqrt(x^2 + y^2 + z^2)
        return ρ_bulge(R,z) + ρ_discs(R,z) + ρ_gas(R,z) + ρ_halo(r)
    end

    return ρ_all(x, y, z)
end

"""
$(TYPEDSIGNATURES)

Compute baryonic Milky Way density at position (x, y, z).
"""
function density_baryon_MW(x, y, z)
    ρ_bulge, ρ_discs, ρ_gas, ρ_halo = GalacticDynamics.milkyway_zhu2023how()

    function ρ_baryon(x, y, z)
        R = sqrt(x^2 + y^2)
        return ρ_bulge(R,z) + ρ_discs(R,z) + ρ_gas(R,z)
    end

    return ρ_baryon(x, y, z)
end

"""
$(TYPEDSIGNATURES)
Generate Milky Way initial conditions.
"""
function generate_milkyway_initial_conditions(xxx, yyy, zzz, r, Δ, unit_cell_volumn, model, FDM_mass_ratio, FDM_radius_ratio, 
                                              baryon_mode, Np, GravitySolver, SofteningLength, 
                                              length_astro, density_astro, potential_astro, acc_astro, 
                                              pids)
    if model == :MW
        model_halo = Zhao(1.55e7u"Msun/kpc^3" * FDM_mass_ratio, 11.75u"kpc" * FDM_radius_ratio, 1.19, 2.95, 0.95)
        ρ_halo = sampling_density.(r * length_astro; model = model_halo) |> collect
        if baryon_mode == :mesh
            ρ_baryon = density_baryon_MW.(xxx*length_astro, yyy*length_astro, zzz*length_astro)/density_astro
            Φ_b = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ρ_baryon, Periodic(), gpu ? GPU() : CPU()))
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
            pos = PVector.(xxx * length_astro, yyy * length_astro, zzz * length_astro)
            particles = generate_milkyway_baryon_particles(Np)
            
            sim_force_baryon = Simulation(particles;
                GravitySolver,
                pids,
            )
            @info "Computing baryonic potentials and forces with $(traitstring(GravitySolver)) solver"
            @time Φ_b = compute_potential(sim_force_baryon, pos, SofteningLength, GravitySolver, CPU()) ./ potential_astro
            @time acc_b = StructArray(compute_force(sim_force_baryon, pos, SofteningLength, GravitySolver, CPU()))
            ax_b = upreferred.(acc_b.x ./ acc_astro)
            ay_b = upreferred.(acc_b.y ./ acc_astro)
            az_b = upreferred.(acc_b.z ./ acc_astro)

            total_mass_baryon = sum(particles.Mass)
            ρ_baryon = nothing
        elseif baryon_mode == :particles_dynamic #TODO
        elseif baryon_mode == :ignored
            ρ_baryon = Φ_b = ax_b = ay_b = az_b = nothing
            total_mass_baryon = 0.0u"Msun"
        end
        
        return ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon
    else
        error("Unsupported model: $(model)")
    end
end

# Internal function for sampling density
function sampling_density(r, model, length_astro, density_astro)
    out_density = upreferred(GalacticDynamics.density(model, r*length_astro) / density_astro)
end

# Generate initial conditions for different models
function generate_initial_conditions(
    model, xxx, yyy, zzz, r, Δ, unit_cell_volumn, FDM_mass_ratio, FDM_radius_ratio,
    baryon_mode, Np, GravitySolver, SofteningLength,
    length_astro, density_astro, potential_astro, acc_astro, pids,
    baryon_β, baryon_ρ0, baryon_r0, halo_β, halo_ρ0, halo_r0, halo_α, halo_γ, halo_Q,
    stellar_TotalMass, stellar_ScaleRadius, thickness_ratio_stellar, gases_TotalMass, gases_ScaleRadius, thickness_ratio_gases
)
    if model == :MW
        ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon = generate_milkyway_initial_conditions(
            xxx, yyy, zzz, r, Δ, unit_cell_volumn, model, FDM_mass_ratio, FDM_radius_ratio, 
            baryon_mode, Np, GravitySolver, SofteningLength, 
            length_astro, density_astro, potential_astro, acc_astro, 
            pids
        )
    elseif model == :SPARC
        #TODO
    elseif model == :cluster_NFW
        model_halo = gNFW(halo_β, halo_ρ0 * FDM_mass_ratio, halo_r0 * FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect
        if baryon_mode == :mesh
            model_baryon = BetaModel(baryon_β, baryon_ρ0, baryon_r0)
            ρ_baryon = sampling_density.(r, model_baryon, length_astro, density_astro) |> collect

            Φ_b = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ρ_baryon, Periodic(), gpu ? GPU() : CPU()))
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
        end
    elseif model == :cluster_Burkert
        model_halo = Burkert(halo_ρ0 * FDM_mass_ratio, halo_r0 * FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect
        if baryon_mode == :mesh
            model_baryon = BetaModel(baryon_β, baryon_ρ0, baryon_r0)
            ρ_baryon = sampling_density.(r, model_baryon, length_astro, density_astro) |> collect
        
            Φ_b = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ρ_baryon, Periodic(), gpu ? GPU() : CPU()))
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
        end
    elseif model == :Elliptical
        model_halo = gNFW(halo_β, halo_ρ0 * FDM_mass_ratio, halo_r0 * FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect
        if baryon_mode == :mesh
            model_baryon = Jaffe(baryon_r0, baryon_ρ0)
            ρ_baryon = sampling_density.(r, model_baryon, length_astro, density_astro) |> collect
        
            Φ_b = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ρ_baryon, Periodic(), gpu ? GPU() : CPU()))
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
        end
    elseif model == :dwarf
        model_halo = gNFW(halo_β, halo_ρ0 * FDM_mass_ratio, halo_r0 * FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect

        if baryon_mode == :ignored
            ρ_baryon = nothing
            total_mass_baryon = 0.0u"Msun"
        elseif baryon_mode == :particles_static
            @info "Sampling stars"
            pos = PVector.(xxx * length_astro, yyy * length_astro, zzz * length_astro)
            
            stellar_ScaleHeight = stellar_ScaleRadius * thickness_ratio_stellar
            config_Stellar = AstroIC.ExponentialDisc(;
                collection = STAR,
                NumSamples = div(Np,2),
                TotalMass = stellar_TotalMass,
                ScaleRadius = stellar_ScaleRadius,
                ScaleHeight = stellar_ScaleHeight,
            )
            particles_Stellar = generate(config_Stellar)

            if isnan(thickness_ratio_gases)
                particles = deepcopy(particles_Stellar)  # Try fix the halting problem
            else
                @info "Sampling gases"
                gases_ScaleHeight = gases_ScaleRadius * thickness_ratio_gases
                config_gases = AstroIC.ExponentialDisc(;
                    collection = STAR,
                    NumSamples = div(Np,2),
                    TotalMass = gases_TotalMass,
                    ScaleRadius = gases_ScaleRadius,
                    ScaleHeight = gases_ScaleHeight,
                )
                particles_gases = generate(config_gases)

                particles = [particles_Stellar; particles_gases]
            end

            sim_force_baryon = Simulation(particles;
                GravitySolver,
                pids,
            )
            @info "Computing baryonic potentials and forces with $(traitstring(GravitySolver)) solver"
            @time Φ_b = compute_potential(sim_force_baryon, pos, SofteningLength, GravitySolver, CPU()) ./ potential_astro
            @time acc_b = StructArray(compute_force(sim_force_baryon, pos, SofteningLength, GravitySolver, CPU()))
            ax_b = upreferred.(acc_b.x ./ acc_astro)
            ay_b = upreferred.(acc_b.y ./ acc_astro)
            az_b = upreferred.(acc_b.z ./ acc_astro)

            total_mass_baryon = sum(particles.Mass)

            ρ_baryon = nothing
            particles_Stellar = particles = nothing
        end
    elseif model == :dwarf_NFW
        model_halo = NFW(halo_ρ0 * FDM_mass_ratio, halo_r0 * FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect

        if baryon_mode == :ignored
            ρ_baryon = nothing
            total_mass_baryon = 0.0u"Msun"
        elseif baryon_mode == :particles_static

        end
    elseif model == :dwarf_Zhao
        model_halo = Zhao(halo_ρ0 * FDM_mass_ratio, halo_r0 * FDM_radius_ratio, halo_α, halo_β, halo_γ)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect

        if baryon_mode == :ignored
            ρ_baryon = nothing
            total_mass_baryon = 0.0u"Msun"
        elseif baryon_mode == :particles_static

        end
    else
        error("Unknown model: $model")
    end

    return ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon
end

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