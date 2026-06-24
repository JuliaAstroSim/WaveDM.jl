# Main IC generation module

# Note: profiles.jl and milkyway.jl are already included at module level in WaveDM.jl

"""
$(TYPEDSIGNATURES)

Generate initial conditions for different models.
"""
function generate_initial_conditions(config_IC::InitialConditionsConfig, grid::SimulationGrid, config_profile::DensityProfileConfig, config_units::AstroUnitsConfig, config_device::DeviceConfig, boundary)
    model = config_IC.model
    r = grid.r
    Δ = grid.Δ
    baryon_mode = config_IC.baryon_mode
    pids = config_IC.pids
    length_astro = config_units.length_astro
    density_astro = config_units.density_astro
    potential_astro = config_units.potential_astro
    acc_astro = config_units.acc_astro
    baryon_β = config_profile.baryon_β
    baryon_ρ0 = config_profile.baryon_ρ0
    baryon_r0 = config_profile.baryon_r0
    halo_β = config_profile.halo_β
    halo_ρ0 = config_profile.halo_ρ0
    halo_r0 = config_profile.halo_r0
    halo_α = config_profile.halo_α
    halo_γ = config_profile.halo_γ
    stellar_TotalMass = config_profile.stellar_TotalMass
    stellar_ScaleRadius = config_profile.stellar_ScaleRadius
    thickness_ratio_stellar = config_profile.thickness_ratio_stellar
    gases_TotalMass = config_profile.gases_TotalMass
    gases_ScaleRadius = config_profile.gases_ScaleRadius
    thickness_ratio_gases = config_profile.thickness_ratio_gases

    ρ_baryon = Φ_b = ax_b = ay_b = az_b = nothing
    total_mass_baryon = 0.0u"Msun"
    baryon_particles = nothing
    if model == :MW
        ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon, baryon_particles = generate_milkyway_initial_conditions(grid, config_IC, config_units, config_device, boundary)
    elseif model == :SPARC_LTGs
        error("Model :SPARC_LTGs is not yet implemented.  Use :MW, :cluster_NFW, :cluster_Burkert, :Elliptical, :dwarf, :dwarf_NFW, or :dwarf_Zhao.")
    elseif model == :SPARC_Xray_ETGs
        error("Model :SPARC_Xray_ETGs is not yet implemented.  Use :MW, :cluster_NFW, :cluster_Burkert, :Elliptical, :dwarf, :dwarf_NFW, or :dwarf_Zhao.")
    elseif model == :SPARC_rotating_ETGs
        error("Model :SPARC_rotating_ETGs is not yet implemented.  Use :MW, :cluster_NFW, :cluster_Burkert, :Elliptical, :dwarf, :dwarf_NFW, or :dwarf_Zhao.")
    elseif model == :cluster_NFW
        model_halo = gNFW(halo_β, halo_ρ0 * config_IC.FDM_mass_ratio, halo_r0 * config_IC.FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect
        if baryon_mode == :mesh
            baryon_particles = nothing
            model_baryon = BetaModel(baryon_β, baryon_ρ0, baryon_r0)
            ρ_baryon = sampling_density.(r, model_baryon, length_astro, density_astro) |> collect

            Φ_b = 4π * parallel_poisson(Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ρ_baryon, boundary, config_device)
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * grid.unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
            error("baryon_mode=:particles_static is not yet implemented for model :cluster_NFW.  Use :mesh or :ignored.")
        else
            error("Unsupported baryon_mode $(baryon_mode) for model :cluster_NFW.  Use :mesh or :ignored.")
        end
    elseif model == :cluster_Burkert
        model_halo = Burkert(halo_ρ0 * config_IC.FDM_mass_ratio, halo_r0 * config_IC.FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect
        if baryon_mode == :mesh
            baryon_particles = nothing
            model_baryon = BetaModel(baryon_β, baryon_ρ0, baryon_r0)
            ρ_baryon = sampling_density.(r, model_baryon, length_astro, density_astro) |> collect

            Φ_b = 4π * parallel_poisson(Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ρ_baryon, boundary, config_device)
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * grid.unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
            error("baryon_mode=:particles_static is not yet implemented for model :cluster_Burkert.  Use :mesh or :ignored.")
        else
            error("Unsupported baryon_mode $(baryon_mode) for model :cluster_Burkert.  Use :mesh or :ignored.")
        end
    elseif model == :Elliptical
        model_halo = gNFW(halo_β, halo_ρ0 * config_IC.FDM_mass_ratio, halo_r0 * config_IC.FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect
        if baryon_mode == :mesh
            baryon_particles = nothing
            model_baryon = Jaffe(baryon_r0, baryon_ρ0)
            ρ_baryon = sampling_density.(r, model_baryon, length_astro, density_astro) |> collect

            Φ_b = 4π * parallel_poisson(Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ρ_baryon, boundary, config_device)
            ax_b, ay_b, az_b = grad_central(-Δ..., Φ_b)
            total_mass_baryon = sum(ρ_baryon) * grid.unit_cell_volumn * density_astro
        elseif baryon_mode == :particles_static
            error("baryon_mode=:particles_static is not yet implemented for model :Elliptical.  Use :mesh or :ignored.")
        else
            error("Unsupported baryon_mode $(baryon_mode) for model :Elliptical.  Use :mesh or :ignored.")
        end
    elseif model == :dwarf
        model_halo = gNFW(halo_β, halo_ρ0 * config_IC.FDM_mass_ratio, halo_r0 * config_IC.FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect

        if baryon_mode == :ignored
            # already set to nothing
            total_mass_baryon = 0.0u"Msun"
        elseif baryon_mode == :particles_static || baryon_mode == :particles_dynamic
            @info "Sampling stars"
            pos = PVector.(grid.xxx * length_astro, grid.yyy * length_astro, grid.zzz * length_astro)

            stellar_ScaleHeight = stellar_ScaleRadius * thickness_ratio_stellar
            config_Stellar = AstroIC.ExponentialDisc(;
                collection = STAR,
                NumSamples = div(config_IC.Np,2),
                TotalMass = stellar_TotalMass,
                ScaleRadius = stellar_ScaleRadius,
                ScaleHeight = stellar_ScaleHeight,
            )
            particles_Stellar = generate(config_Stellar)

            if isnan(thickness_ratio_gases)
                baryon_particles = deepcopy(particles_Stellar)  # Try fix the halting problem
            else
                @info "Sampling gases"
                gases_ScaleHeight = gases_ScaleRadius * thickness_ratio_gases
                config_gases = AstroIC.ExponentialDisc(;
                    collection = STAR,
                    NumSamples = div(config_IC.Np,2),
                    TotalMass = gases_TotalMass,
                    ScaleRadius = gases_ScaleRadius,
                    ScaleHeight = gases_ScaleHeight,
                )
                particles_gases = generate(config_gases)

                baryon_particles = [particles_Stellar; particles_gases]
            end

            sim_force_baryon = Simulation(baryon_particles;
                GravitySolver = config_IC.GravitySolver,
                pids,
            )
            @info "Computing baryonic potentials and forces with $(traitstring(config_IC.GravitySolver)) solver"
            @time Φ_b = compute_potential(sim_force_baryon, pos, config_IC.SofteningLength, config_IC.GravitySolver, CPU()) ./ potential_astro
            @time acc_b = StructArray(compute_force(sim_force_baryon, pos, config_IC.SofteningLength, config_IC.GravitySolver, CPU()))
            ax_b = upreferred.(acc_b.x ./ acc_astro)
            ay_b = upreferred.(acc_b.y ./ acc_astro)
            az_b = upreferred.(acc_b.z ./ acc_astro)

            total_mass_baryon = sum(baryon_particles.Mass)

            ρ_baryon = nothing
            particles_Stellar = nothing
        else
            error("Unsupported baryon_mode $(baryon_mode) for model :dwarf.  Use :ignored, :particles_static, or :particles_dynamic.")
        end
    elseif model == :dwarf_NFW
        model_halo = NFW(halo_ρ0 * config_IC.FDM_mass_ratio, halo_r0 * config_IC.FDM_radius_ratio)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect

        if baryon_mode == :ignored
            # already set to nothing
            total_mass_baryon = 0.0u"Msun"
        elseif baryon_mode == :particles_static
            error("baryon_mode=:particles_static is not yet implemented for model :dwarf_NFW.  Use :ignored.")
        else
            error("Unsupported baryon_mode $(baryon_mode) for model :dwarf_NFW.  Use :ignored.")
        end
    elseif model == :dwarf_Zhao
        model_halo = Zhao(halo_ρ0 * config_IC.FDM_mass_ratio, halo_r0 * config_IC.FDM_radius_ratio, halo_α, halo_β, halo_γ)
        ρ_halo = sampling_density.(r, model_halo, length_astro, density_astro) |> collect

        if baryon_mode == :ignored
            # already set to nothing
            total_mass_baryon = 0.0u"Msun"
        elseif baryon_mode == :particles_static
            error("baryon_mode=:particles_static is not yet implemented for model :dwarf_Zhao.  Use :ignored.")
        else
            error("Unsupported baryon_mode $(baryon_mode) for model :dwarf_Zhao.  Use :ignored.")
        end
    else
        error("Unknown model: $model")
    end

    @info "Computing WaveDM acc"
    @time Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM = compute_acceleration_field(ρ_halo, grid, boundary, config_device.gpu, config_units.mass_astro, config_IC.SofteningLength, potential_astro, pids, length_astro)

    @info "Setting WaveDM vel"
    v, Φ_all, ax_all, ay_all, az_all = generate_velocity_field(Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM, Φ_b, ax_b, ay_b, az_b, baryon_mode, grid, config_IC)
    vx, vy, vz = adjust_velocity_field(v, ρ_halo, config_IC.bulk_perturb, config_IC.bulk_size, config_IC.bulk_shift_size, config_IC.bulk_center_size)

    θ = solve_vector_equation(collect(vx), collect(vy), collect(vz), Δ...)
    phase = exp.(im*θ)
    IC_vel = sqrt.(ρ_halo) .* phase

    if baryon_mode == :ignored
    elseif baryon_mode == :particles_static
    elseif baryon_mode == :particles_dynamic
        config_mesh = MeshConfig(uAstro;
            mode = VertexMode(),
            assignment = CIC(),
            # assignment = TSC(),
            # boundary = Periodic(),
            Nx = grid.Nx,
            Ny = grid.Ny,
            Nz = grid.Nz,
            NG = 0,
            xMin = -0.5*grid.Xmax * length_astro,
            xMax = +0.5*grid.Xmax * length_astro,
            yMin = -0.5*grid.Ymax * length_astro,
            yMax = +0.5*grid.Ymax * length_astro,
            zMin = -0.5*grid.Zmax * length_astro,
            zMax = +0.5*grid.Zmax * length_astro,
            dim = 3,
            device = CPU(),
        )
        meshpos = PVector.(grid.xxx * length_astro, grid.yyy * length_astro, grid.zzz * length_astro)
        meshacc = PVector.(ax_all * acc_astro, ay_all * acc_astro, az_all * acc_astro)

        if model == :MW && config_IC.MW_disk_RC
            @info "Setting baryon vel: bulge"
            for i in eachindex(baryon_particles)
                if baryon_particles[i].collection == BULGE
                    pos = baryon_particles[i].Pos
                    if !is_inbound(pos, config_mesh)
                        continue
                    end
                    acc = uconvert(u"kpc/Gyr^2", mesh2particle(meshpos, config_mesh, meshacc, pos, config_mesh.mode, config_mesh.assignment))
                    vel = uconvert(u"kpc/Gyr", AstroIC.rotational_velocity_acc(pos.x, pos.y, pos.z, norm(acc), config_IC.rotational_ratio_baryon) * config_IC.velocity_ratio_baryon)
                    setproperty!!(baryon_particles[i], :Acc, acc)
                    setproperty!!(baryon_particles[i], :Vel, vel)
                end
            end
        else
            @info "Setting baryon vel: all"
            for i in eachindex(baryon_particles)
                pos = baryon_particles[i].Pos
                if !is_inbound(pos, config_mesh)
                    continue
                end
                acc = uconvert(u"kpc/Gyr^2", mesh2particle(meshpos, config_mesh, meshacc, pos, config_mesh.mode, config_mesh.assignment))
                vel = uconvert(u"kpc/Gyr", AstroIC.rotational_velocity_acc(pos.x, pos.y, pos.z, norm(acc), config_IC.rotational_ratio_baryon) * config_IC.velocity_ratio_baryon)
                setproperty!!(baryon_particles[i], :Acc, acc)
                setproperty!!(baryon_particles[i], :Vel, vel)
            end
        end
    end

    return IC_vel, ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon, baryon_particles
end

# Compute acceleration field from potential
function compute_acceleration_field(ρ_halo, grid::SimulationGrid, boundary, gpu, mass_astro, SofteningLength, potential_astro, pids, length_astro)
    Nx, Ny, Nz = size(ρ_halo)
    if boundary isa Vacuum
        # Use octree to compute gravitational force
        mesh_particles = StructArray(Star(uAstro; id = i+(j-1)*Nx+(k-1)*Nx*Ny) for i in 1:Nx, j in 1:Ny, k in 1:Nz)
        mesh_particles.Pos .= PVector.(grid.xxx * length_astro, grid.yyy * length_astro, grid.zzz * length_astro)
        mesh_particles.Mass .= ρ_halo * grid.unit_cell_volumn * mass_astro

        sim_mesh_force = Simulation(mesh_particles;
            GravitySolver = Tree(),
            pids,
        )

        Φ_WaveDM = compute_potential(sim_mesh_force, mesh_particles.Pos, SofteningLength, Tree(), CPU()) ./ potential_astro
    else
        # No boundary-specific `config_device` is in scope here; build a
        # throwaway backend from the `gpu` hint.  In the WaveDM main
        # loop this code path is replaced by the typed overload
        # `parallel_poisson(..., config_device)` which sees the full backend state.
        cfg = DeviceConfig(; gpu=gpu)
        Φ_WaveDM = 4π * parallel_poisson(grid.Δ, [Nx-1, Ny-1, Nz-1], ρ_halo, Periodic(), cfg)
    end

    ax_WaveDM, ay_WaveDM, az_WaveDM = grad_central(-grid.Δ..., Φ_WaveDM)

    return Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM
end

# Generate velocity field based on acceleration
function generate_velocity_field(Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM, Φ_b, ax_b, ay_b, az_b, baryon_mode, grid::SimulationGrid, config_IC::InitialConditionsConfig)
    Nx, Ny, Nz = size(grid.xxx)
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

    if config_IC.velocity_falling
        v = AstroIC.freefall_velocity_acc.(grid.xxx, grid.yyy, grid.zzz, a_CDM) .* config_IC.velocity_ratio |> collect
    else
        v = AstroIC.rotational_velocity_acc.(grid.xxx, grid.yyy, grid.zzz, a_CDM, config_IC.rotational_ratio) .* config_IC.velocity_ratio |> collect
    end

    return v, Φ_all, ax_all, ay_all, az_all
end

# Adjust velocity field (center velocity, bulk perturbation, etc.)
function adjust_velocity_field(v, ρ_halo, bulk_perturb, bulk_size, bulk_shift_size, bulk_center_size)
    Nx, Ny, Nz = size(v)
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

function optimize_smoothed_inversion(r, massBaryon, lower, upper, ic;
    mode = :NFW, # support :NFW, TODO :nonparameteric
)
    function target_smoothed_inversion(params)
        if mode == :NFW
            rho0 = params[1] * u"Msun/kpc^3" * 1e8
            r_s = params[2] * u"kpc"
            modelDM = NFW(rho0, r_s)
            rhoDM = GalacticDynamics.density.(modelDM, r)
        else
            error("Unsupported DM mode! Try keyword `mode=:NFW`")
        end
        massDM = 4π * cumul_integrate(r, r.^2 .* rhoDM)
        massTotal = massBaryon + massDM
        constraint_nonmonotonic = sum(ustrip.(diff(massTotal)) .< 0) # This is loose constraint and not sufficient to optimize the DM halo

        accTotal = ustrip.(u"m/s^2", C.G .* massTotal ./ r.^2)
        accBaryon = ustrip.(u"m/s^2", C.G .* massBaryon ./ r.^2)
        accRAR = RAR.(accBaryon, 1.2e-10)
        constraint_acc = sum(((accTotal[div(end,2):end] .- accRAR[div(end,2):end])*1e11).^2) / length(accRAR) # mean squared error, MSE. Enlarge the unit for more strict constraint
        # Should not use the whole array

        #TODO: add weights
        return constraint_nonmonotonic + constraint_acc
    end

    result = Optim.optimize(
        target_smoothed_inversion,
        lower, upper, ic,
        Optim.Fminbox(),
        Optim.Options(
            store_trace = true,
            # iterations = 100,
            # outer_iterations = 100,
            # x_tol = 1e-10,
            # outer_x_tol = 1e-10,
        )
    )
    # trace = [state.value for state in result.trace]
    # @show trace
    # println(UnicodePlots.lineplot(trace))

    minimum_value = Optim.minimum(result)
    optimal_params = Optim.minimizer(result)
    @show minimum_value, optimal_params

    return optimal_params
end

function optimize_smoothed_inversion_forward_fitting(r, massBaryon, r_RC, v_RC, lower, upper, ic;
    mode = :NFW, # support :NFW, TODO :nonparameteric
)
    function target_smoothed_inversion_forward_fitting(params)
        if mode == :NFW
            rho0 = params[1] * u"Msun/kpc^3"
            r_s = params[2] * u"kpc"
            modelDM = NFW(rho0, r_s)
            rhoDM = GalacticDynamics.density.(modelDM, r)
        else
            error("Unsupported DM mode! Try keyword `mode=:NFW`")
        end
        massDM = 4π * cumul_integrate(r, r.^2 .* rhoDM)
        massTotal = massBaryon + massDM

        constraint_nonmonotonic = sum(diff(massTotal) .< zero(eltype(massTotal)))

        velTotal = ustrip.(u"km/s", sqrt.(C.G .* massTotal ./ r))
        spl_velTotal = Spline1D(ustrip.(u"kpc", r), velTotal)
        constraint_RC = sum((spl_velTotal(ustrip.(u"kpc", r_RC)) .- ustrip.(u"km/s", v_RC)) .^ 2) / length(r_RC)

        accTotal = ustrip.(u"m/s^2", C.G .* massTotal ./ r.^2)
        accBaryon = ustrip.(u"m/s^2", C.G .* massBaryon ./ r.^2)
        accRAR = RAR.(accBaryon, 1.2e-10)
        constraint_acc = sum(((accTotal[div(end,2):end] .- accRAR[div(end,2):end])*1e11).^2) / length(accRAR) # mean squared error, MSE. Enlarge the unit for more strict constraint

        #TODO: add weights
        # return constraint_nonmonotonic + 1 * constraint_RC + constraint_acc
        return constraint_nonmonotonic + constraint_acc
        # return constraint_acc
    end
    
    result = Optim.optimize(
        target_smoothed_inversion_forward_fitting,
        lower, upper, ic,
        Optim.Fminbox(),
        Optim.Options(
            store_trace = true,
            iterations = 500,
            outer_iterations = 500,
            # x_tol = 1e-10,
            # outer_x_tol = 1e-10,
        )
    )
    # trace = [state.value for state in result.trace]
    # println(UnicodePlots.lineplot(trace))
    # @show trace

    minimum_value = Optim.minimum(result)
    optimal_params = Optim.minimizer(result)
    @show minimum_value, optimal_params

    return optimal_params
end