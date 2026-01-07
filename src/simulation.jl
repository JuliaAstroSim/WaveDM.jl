
function SPE3D_MOND(;
    Xmax = 5,
    Ymax = Xmax,
    Zmax = Xmax,
    Tmax = 1,
    Nx = 64,
    Ny = Nx,
    Nz = Nx,
    Nt = 32,
    autoset_timestep = false,
    autoset_timestep_ratio = 0.9,
    IC = (x,y,z)->exp(-x^2-y^2-z^2), # Initial conditions of wave function. Provide a matrix or distribution function f(x,y,z). Default: Gaussian profile.
    V = (x,y,z,ψ)->0.0f0, # Potentials. Provide a function f(x,y,z,ψ). Default: free propagation
    κ = 0,
    baryon = (x,y,z,ψ)->0, # baryon term 4πGρ
    baryon_mode = :mesh, # Options: `:mesh`, sample baryons on the mesh; `:particles_static`, sample baryons as N-body particles; `ignored`, ignore baryonic effects
    baryon_potential = nothing, # pre-computed potentials of static baryons
    ax_b = nothing,
    ay_b = nothing,
    az_b = nothing,
    absorb_coeff = 0,
    size = (2400, 1400),
    outputdir = joinpath(@__DIR__, "output"),
    title = "MOND3D",
    filename = title,
    boundary = Periodic(),
    # SofteningLength_ratio = 0.5, # relative to Δx
    SofteningLength = 1.0u"kpc",
    gpu = true,
    plotOptical = false,
    plotWaveDM = false,
    Realtime = true,
    StepsBetweenSnapshots = 5,
    FDM_mass_ratio = 1.0,
    FDM_radius_ratio = 1.0,
    rotational_ratio = 0.0,
    minR = 8.0u"kpc",
    maxR = 15.0u"kpc",
    massRadius = 70u"kpc",
    massMW = 4.2e11u"Msun",
    mesh_step::Int = 2,
    zoom_a_max = 2e-10,
    cylindrical = true,
    acc_RAR_min = 1e-11,
    acc_RAR_max = 1e-8,
    nuIndex = 2.0,
    flagRC = true,
    save_phi = true,
    save_IC = true,
    
    best_fit_halo_mass = false,
    KDK_flag = true,
    # dynamic_colorrange = false,
    dynamic_colorrange = true,
    GravitySolver = Tree(),
    unicode_plot = false,
    unicode_heatmap_width = 150,

    mₐ = 0.2 * 1e-22 * 1.783e-36u"kg",
    aₛ = +8.29e-60u"fm", # FDM,
    Ωₘ₀ = 0.31,

    length_astro = uconvert(u"kpc", (8 * π * (C.h/2/π)^2 / (3 * mₐ^2 * C.H^2 * Ωₘ₀))^0.25),
    time_astro = uconvert(u"Gyr", (3 * C.H^2 * Ωₘ₀ / (8 * π))^-0.5),
    mass_astro = uconvert(u"Msun", (3 * C.H^2 * Ωₘ₀ / (8 * π))^0.25 * (C.h/2/π)^1.5 / (mₐ^1.5 * C.G)),
    density_astro = uconvert(u"Msun/kpc^3", mass_astro/length_astro^3),
    acc_astro = uconvert(u"m/s^2", length_astro / time_astro^2),
    velocity_astro = uconvert(u"m/s", length_astro / time_astro),
    potential_astro = uconvert(u"kpc^2/Gyr^2", C.G * mass_astro / length_astro),

    h_astro = C.h / length_astro^2 / mass_astro * time_astro |> upreferred,
    aₛ_astro = aₛ / length_astro |> upreferred,
    mₐ_astro = mₐ / mass_astro |> upreferred,
    c_astro = C.c / length_astro * time_astro |> upreferred,
    κ_astro = 4π * ħ * aₛ / (time_astro * mₐ^2 * C.G) |> upreferred,

    G0 = C.G * mass_astro * time_astro^2 / length_astro^3,
    a0 = upreferred(C.ACC0 / acc_astro),

    plot_virial = false,

    uniform_interval = true,

    ## ultra faint dwarfs: 
    extract_dwarf_granule = false,
    extract_min_t = 0.2u"Gyr",
    target_velocity_dispersion = [],
    target_velocity_dispersion_r = [],
    target_profile_model = :dwarf_Zhao,  # default: Segue I
    target_profile_ρ0 = 10^(-0.9) * u"Msun/pc^3",
    target_profile_ρ0_u = 10^(-0.9+2.4) * u"Msun/pc^3",
    target_profile_ρ0_d = 10^(-0.9-2.5) * u"Msun/pc^3",
    target_profile_rs = 10^(2.5) * u"pc",
    target_profile_rs_u = 10^(2.5+1.1) * u"pc",
    target_profile_rs_d = 10^(2.5-1.4) * u"pc",
    target_profile_α = 1.6,
    target_profile_α_u = 1.6+0.8,
    target_profile_α_d = 1.6-0.9,
    target_profile_β = 6.2,
    target_profile_β_u = 6.2+2.1,
    target_profile_β_d = 6.2-2.3,
    target_profile_γ = 1.6,
    target_profile_γ_u = 1.6+0.6,
    target_profile_γ_d = 1.6-0.3,
    target_profile_Q = 1.0,
    target_profile_error = false,
    target_fitting_rs_ratio = 2,

    extract_mode = :profile,
    folder_data = @__DIR__,

    # average_N = 10, # average on the last few timesteps. If zero, simply take the last timestep for later computation
    average = false, # Reynolds time-averaging
    average_start_t = extract_min_t,
    Galaxy_i = 0,
    average_snapshots = false,

    average_all = false, # spherical shell
    average_all_start_t = average_start_t,
    filter_half_width = 0,
    relative_rho_min_threshold = 0.0,
    target_binning_rs_ratio = 4,

    target_beta_star = NaN,
    target_beta_star_u = NaN, #TODO
    target_beta_star_d = NaN,
    target_beta_star_r_min = 0.3u"kpc",
    target_beta_star_r_max = 0.8u"kpc",
    beta_star_error_threshold = 0.1,

    MW_tidal_field = false,
    tidal_BG_sim = nothing,
    tidal_initial_pos = PVector(-19.4, -9.5, 17.7, u"kpc"), # Galactocentric
    tidal_initial_vel = PVector(13.0, -175.0, 51.0, u"km/s"),
    tidal_lookback_time = 0.0u"Gyr",
    df_traj = DataFrame(),
    spl_pot = nothing,
    sim_force_baryon = nothing, # MW
    particles_LMC = nothing,
    df_traj_LMC = nothing,

    LMC_tidal_field = false,

    MW_tidal_interpolate = true,
    MW_pot = nothing, # if is a filename string, load the variable MW_x, MW_y, MW_z, MW_Phi with astro units; if nothing, generate a new one; if 3D Array, use it immediately
    MW_pot_Xmax = 100.0u"kpc",
    MW_pot_Ymax = 100.0u"kpc",
    MW_pot_Zmax = 100.0u"kpc",
    MW_pot_N = 512,
    export_MW_pot = false,

    # distributed_memory = true,
    distributed_memory = false,
    pids = workers(),
    kw...
)
    println("\n\n")
    mkpathIfNotExist(outputdir)
    
    if distributed_memory
        @info "Distributed memory parallelism enabled!"
        DA = DArray
    else
        DA = collect
    end

    @info gpu ? "Carrying out FFT on GPU" : "Carrying out FFT on CPU"
    DeviceArray = gpu ? cu : collect # convertor

    @info "TimeMax: $(Tmax * time_astro)"
    @info "Xmax: $(Xmax * length_astro)"

    uT = ustrip(u"Gyr", time_astro)
    uL = ustrip(u"kpc", length_astro)
    uVel = ustrip(u"km/s", velocity_astro)
    uAcc = ustrip(u"m/s^2", acc_astro)
    uRho = ustrip(u"Msun/kpc^3", density_astro)
    uMomentum = ustrip(u"Msun*kpc/Gyr", mass_astro * velocity_astro)

    x, y, z, Δ, unit_cell_volumn = setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)
    section = ceil(Int, Nx/2*sqrt(3))
    oneMatrix = ones(Nx, Ny, Nz)
    xxx, yyy, zzz, r = setup_coordinates(x, y, z, Nx, Ny, Nz, oneMatrix; DA)

    @info "Sim v_max = $(uconvert(u"km/s", π * C.h/(2π) / (mₐ_astro * mass_astro) / (Δ[1]*length_astro)))"
    @info "Setting softening length to $(SofteningLength)"
    @info "Setting ICs"
    ψ, sqrt_rho, rho = setup_initial_conditions(IC, xxx, yyy, zzz; DA = collect)

    if baryon_mode == :mesh # sample baryons on the mesh
        if isnothing(baryon_potential)
            if baryon isa Function
                baryon_term = baryon.(xxx, yyy, zzz)
            else
                baryon_term = baryon
            end
            Φ_b = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], DeviceArray(baryon_term), Periodic(), gpu ? GPU() : CPU()))
        else
            Φ_b = baryon_potential
        end
    elseif baryon_mode == :particles_static # sample baryons as N-body particles
        if isnothing(baryon_potential)
            # try solve the potential from particles
            # pos = PVector.(xxx * length_astro, yyy * length_astro, zzz * length_astro)
            # Φ_b = [AstroNbodySim.compute_unit_potential_at_point(p, particles, C.G, 0.01u"kpc") for p in pos]) / potential_astr
            @error "baryon_potential should be provided"
        else
            Φ_b = baryon_potential
        end
    elseif baryon_mode == :particles_dynamic #TODO
    elseif baryon_mode == :ignored
        Φ_b = nothing
    end

    @info "Solving IC potential"
    sim_mesh_force = nothing
    mesh_particles = nothing
    if boundary isa Vacuum
        # Use octree to compute gravitational force
        mesh_particles = StructArray(Star(uAstro; id = i+(j-1)*Nx+(k-1)*Nx*Ny) for i in 1:Nx, j in 1:Ny, k in 1:Nz)
        mesh_particles.Pos .= PVector.(xxx * length_astro, yyy * length_astro, zzz * length_astro)
        mesh_particles.Mass .= rho * unit_cell_volumn * mass_astro
        
        sim_mesh_force = Simulation(mesh_particles;
            GravitySolver,
            # pids, #TODO parallelism need bcast particle mass
        )
        # mesh_particles = nothing # release memory
        Φ_WaveDM = compute_potential(sim_mesh_force, mesh_particles.Pos, SofteningLength, Tree(), CPU()) ./ potential_astro
    else
        Φ_WaveDM = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], abs.(DeviceArray(ψ)).^2, Periodic(), gpu ? GPU() : CPU()))
    end

    @info "Computing IC acc"
    ax_WaveDM, ay_WaveDM, az_WaveDM = grad_central(-Δ..., Φ_WaveDM)

    if baryon_mode == :ignored
        Φ_all = Φ_WaveDM
        ax_all = ax_WaveDM
        ay_all = ay_WaveDM
        az_all = az_WaveDM
        a_all = sqrt.(ax_all[:, :, div(end,2)].^2 .+ ay_all[:, :, div(end,2)].^2 .+ az_all[:, :, div(end,2)].^2)
        ax_all = ay_all = az_all = nothing
    else
        Φ_all = Φ_WaveDM + baryon_potential
        ax_all = ax_WaveDM + ax_b
        ay_all = ay_WaveDM + ay_b
        az_all = az_WaveDM + az_b
        ax_WaveDM = ay_WaveDM = az_WaveDM = nothing
        a_all = sqrt.(ax_all[:, :, div(end,2)].^2 .+ ay_all[:, :, div(end,2)].^2 .+ az_all[:, :, div(end,2)].^2)
    end

    # z=0 plane
    if baryon_mode != :ignored
        a_b = sqrt.(ax_b[:, :, div(end,2)].^2 .+ ay_b[:, :, div(end,2)].^2 .+ az_b[:, :, div(end,2)].^2)
    end

    t, dt, Nt = compute_timestep(Δ[1], Φ_all, κ, ψ, Tmax, Nt, autoset_timestep, autoset_timestep_ratio)
    @info "Default Δt: $(Tmax / Nt)"
    @info "Current Δt: $(Tmax / Nt)"

    suffix = "Nx($(@sprintf("%d", Nx))), Xmax($(@sprintf("%.2f", Xmax))), Nt($(@sprintf("%d", Nt))), Tmax($(@sprintf("%.2f", Tmax))), DM_m($(@sprintf("%.2f", FDM_mass_ratio)))"
    @info "Initializing simulation: $(title), $(suffix)"
    @info "Data saved to folder: $(outputdir)"

    @info "Plotting IC"
    if baryon_mode != :ignored
        plotMOND(ax_all, ay_all, az_all, ax_b, ay_b, az_b, a0, r, length_astro, acc_astro, minR, maxR, outputdir, title, suffix, section; filename = "$(title)_IC")
    end

    # rho_max, rho_max_id = findmax(rho) # Too slow for Dagger
    # rho_max, rho_max_id = findmax(collect(rho))
    sum_rho = sum(rho)
    xc = sum(xxx .* rho) / sum_rho
    yc = sum(yyy .* rho) / sum_rho
    zc = sum(zzz .* rho) / sum_rho
    id_xc = findfirstvalue(x, xc)
    id_yc = findfirstvalue(y, yc)
    id_zc = findfirstvalue(z, zc)
    rho_max_id = CartesianIndex(id_xc, id_yc, id_zc)
    rho_max = rho[rho_max_id]

    r_mass_center = sqrt.((xxx.-xc).^2 + (yyy.-yc).^2 + (zzz.-zc).^2) |> collect
    r_filter = r_mass_center .<= target_profile_rs * target_fitting_rs_ratio / length_astro
    
    total_halo_mass = sum(rho) * unit_cell_volumn
    radii = quantile(r_mass_center[:], weights(rho[:]), 0.1:0.1:0.9)

    @info("Start visualization...")
    astro_config = AstroUnitsConfig(length_astro, time_astro, mass_astro, density_astro, acc_astro, velocity_astro, potential_astro, 
                                    uT, uL, uVel, uAcc, uRho, uMomentum, h_astro, aₛ_astro, mₐ_astro, c_astro, κ_astro, G0, a0)
    grid = SimulationGrid(Xmax, Ymax, Zmax, Nx, Ny, Nz, Δ, x, y, z, xxx, yyy, zzz, r, oneMatrix, unit_cell_volumn)
    vis_config = VisualizationConfig(title, suffix, size, StepsBetweenSnapshots, Realtime, dynamic_colorrange, plot_virial, plotOptical, plotWaveDM)
    data_config = VisualizationData(rho, rho_max_id, total_halo_mass, radii, r_mass_center, target_profile_model, target_profile_error)
    
    fig, ArrayT, ArrayT_Snap, AxisR, AxisVirial, AxisDensityProfile, SliceXY, SliceYZ, SliceXZ, ArrayTotalMass, ArrayR, ArrayR1, ArrayR2, ArrayR3, ArrayR4, ArrayR5, ArrayR6, ArrayR7, ArrayR8, ArrayR9, ColorRange = setup_visualization(
        grid, t, vis_config, data_config, astro_config, distributed_memory
    )

    if baryon_mode == :ignored
    else
        a_b = sqrt.(ax_b[:, :, div(end,2)].^2 .+ ay_b[:, :, div(end,2)].^2 .+ az_b[:, :, div(end,2)].^2)
        a_all = sqrt.(ax_all[:, :, div(end,2)].^2 .+ ay_all[:, :, div(end,2)].^2 .+ az_all[:, :, div(end,2)].^2)
        a_mond = a_b ./ (1 .- exp.(-sqrt.(a_b./a0))) #RAR
        rMOND = r[:, :, div(end,2)][:] * ustrip(length_astro)
    end

    if plot_virial
        @info "Computing momentums and virial energies"
        ArrayVirialPotential, ArrayTotalKineticE, ArrayTotalQuantumE, ArrayVirial, ArrayMomentumX, ArrayMomentumY, ArrayMomentumZ = setup_virial_visualization(
            ψ, Φ_all, rho, sqrt_rho, Δ, unit_cell_volumn, mass_astro, velocity_astro, ArrayT, AxisVirial
        )
    else
        ArrayVirialPotential = ArrayTotalKineticE = ArrayTotalQuantumE = ArrayVirial = ArrayMomentumX = ArrayMomentumY = ArrayMomentumZ = nothing
    end

    profile_r_mean, profile_ρ_mean, r_target, ρ_halo_target = setup_density_profile_visualization!(
        fig, AxisDensityProfile, Δ, rho, r_filter, r_mass_center, target_profile_model, target_profile_error, target_profile_ρ0, target_profile_ρ0_u, target_profile_ρ0_d, target_profile_rs, target_profile_rs_u, target_profile_rs_d, target_profile_α, target_profile_α_u, target_profile_α_d, target_profile_β, target_profile_β_u, target_profile_β_d, target_profile_γ, target_profile_γ_u, target_profile_γ_d, target_fitting_rs_ratio, length_astro, uL, uRho
    )

    Makie.save(joinpath(outputdir, "$(title), $(suffix) - Overview IC.png"), fig)

    save_IC && save_initial_conditions(ψ, baryon_mode, baryon_potential, ax_b, ay_b, az_b, outputdir, title, suffix)
        
    if Realtime
        display(fig)
    end

    # initialize array
    ArrayTotalMass_temp = Float64[]
    ArrayR_temp = Float64[]
    ArrayR1_temp = Float64[]
    ArrayR2_temp = Float64[]
    ArrayR3_temp = Float64[]
    ArrayR4_temp = Float64[]
    ArrayR5_temp = Float64[]
    ArrayR6_temp = Float64[]
    ArrayR7_temp = Float64[]
    ArrayR8_temp = Float64[]
    ArrayR9_temp = Float64[]
    ArrayVirialPotential_temp = Float64[]
    ArrayTotalKineticE_temp = Float64[]
    ArrayTotalQuantumE_temp = Float64[]
    ArrayVirial_temp = Float64[]

    ArrayMomentumX_temp = Float64[]
    ArrayMomentumY_temp = Float64[]
    ArrayMomentumZ_temp = Float64[]

    linear_phase = setup_fft_operators(Xmax, Ymax, Zmax, Nx, Ny, Nz, dt)  # Laplacian in Fourier space
    boarder = setup_absorption_boundary(Xmax, Ymax, Zmax, x, y, z, absorb_coeff, dt)

    device_ψ = DeviceArray(ψ)

    average_N = 0
    if average
        @info "Averaging ψ from $(average_start_t)"
        buffer_ψ2 = zeros(ComplexF64, Nx, Ny, Nz) # save the last few steps for averaging
    else
        buffer_ψ2 = nothing
    end

    best_fit_error = Inf
    current_fit_error = 0
    best_fit_beta_star_error = Inf
    best_fit_beta_star = NaN
    current_beta_star = NaN
    best_fit_t = 0.0u"Gyr"
    
    if extract_dwarf_granule
        @info "Extracting quasi condensates"
        best_fit_ψ = similar(ψ)
        best_fit_ψ_last_t = similar(ψ)
        best_fit_Φ_all = similar(Φ_all)

        best_fit_time_file = open(joinpath(outputdir, "$(title), $(suffix) - best_fit_time.csv"), "a")

        if !isnan(target_beta_star)
            @info "Finding the best β* within ($(target_beta_star_r_min), $(target_beta_star_r_max))"
        end

        if extract_mode == :RC
            @info "extracting snapshot that best-fitting the CO rotation curves"
            if Galaxy_i < 1 || Galaxy_i > 6
                error("Keyword Galaxy_i has to be integers ∈[1,6]")
            end
            df_CO_RC = DataFrame(CSV.File(joinpath(folder_data, "Cooke2022_RC", "$(dfDwarf.Galaxy[Galaxy_i]).csv")))
            df_CO_RC[!,"vel_e"] = df_CO_RC.vel_u .- df_CO_RC.vel
            df_CO_RC[!,"vel_d"] = df_CO_RC.vel .- df_CO_RC.vel_e
        elseif extract_mode == :profile
            @info "extracting snapshot that best-fitting the density profile"
        end

        profile_config = ProfileFitConfig(target_profile_ρ0, target_profile_ρ0_u, target_profile_ρ0_d, target_profile_rs, target_profile_rs_u, target_profile_rs_d, 
                                                target_profile_α, target_profile_α_u, target_profile_α_d, target_profile_β, target_profile_β_u, target_profile_β_d, 
                                                target_profile_γ, target_profile_γ_u, target_profile_γ_d, target_fitting_rs_ratio, uniform_interval)
        rc_config = RCFitConfig(target_beta_star, target_beta_star_u, target_beta_star_d, target_beta_star_r_min, target_beta_star_r_max, beta_star_error_threshold, Galaxy_i)
    else
        best_fit_t = nothing
    end

    MW_grid = MW_Phi = nothing
    sim_traj_LMC = nothing
    if MW_tidal_field
        if LMC_tidal_field
            @info "Considering time-dependent tidal forces from LMC"
            sim_traj_LMC = Simulation(deepcopy(StructArray(particles_LMC)); GravitySolver, pids);
        end

        if MW_tidal_interpolate
            MW_grid, MW_Phi, MW_x, MW_y, MW_z = setup_mw_tidal_field(MW_pot, MW_pot_Xmax, MW_pot_Ymax, MW_pot_Zmax, MW_pot_N, length_astro, potential_astro, spl_pot, sim_force_baryon, SofteningLength)
        else # Directly compute the potentials
        end
    end
    
    rho_max, rho_max_id = findmax(rho)

    @info "Starting main loop"
    progress = Progress(Nt-1)
    breakflag = false
    ψ_last_t = similar(ψ)
    
    grid = SimulationGrid(Xmax, Ymax, Zmax, Nx, Ny, Nz, Δ, x, y, z, xxx, yyy, zzz, r, oneMatrix, unit_cell_volumn)
    gravity_config = GravityConfig(boundary, Φ_b, sim_mesh_force, mesh_particles, SofteningLength, baryon_mode, GravitySolver, mass_astro)
    device_config = DeviceConfig(gpu, DeviceArray, DA)
    tidal_config = TidalFieldConfig(MW_tidal_field, MW_tidal_interpolate, LMC_tidal_field, uT, tidal_lookback_time, df_traj, df_traj_LMC, length_astro, uL, MW_grid, MW_Phi, spl_pot, sim_force_baryon, particles_LMC, sim_traj_LMC)
    
    Makie.record(fig, joinpath(outputdir, "$(title), $(suffix) - Overview.mp4")) do io
        for i in 2:Nt
            ### Kick
            dt_kick = KDK_flag ? 0.5 * dt : dt
            
            # Pass individual parameters to apply_kick_step! instead of KickStepConfig
            Φ_all, spec = apply_kick_step!(device_ψ, ψ, V, rho_max, rho_max_id, grid, gravity_config, tidal_config, device_config, dt_kick, i, t)
            
            ### Drift
            device_ψ = apply_drift_step!(spec, linear_phase, boarder, gpu, DeviceArray)

            ### Update data handlers
            ψ_last_t .= ψ
            ψ = collect(device_ψ) |> DA  # update
            sqrt_rho = abs.(ψ)
            rho = sqrt_rho.^2
            rho_max, rho_max_id = findmax(rho)
            r_mass_center = sqrt.((xxx.-xxx[rho_max_id]).^2 + (yyy.-yyy[rho_max_id]).^2 + (zzz.-zzz[rho_max_id]).^2) |> collect

            push!(ArrayTotalMass_temp, isinf(total_halo_mass) ? 0 : total_halo_mass)
            push!(ArrayR_temp, std(collect(r[:,:,div(end,2)]), aweights(collect(abs.(ψ[:,:,div(end,2)]).^2))) * uL)

            radii = quantile(r_mass_center[:], weights(rho[:]), 0.1:0.1:0.9)
            push!(ArrayR1_temp, radii[1] * uL)
            push!(ArrayR2_temp, radii[2] * uL)
            push!(ArrayR3_temp, radii[3] * uL)
            push!(ArrayR4_temp, radii[4] * uL)
            push!(ArrayR5_temp, radii[5] * uL)
            push!(ArrayR6_temp, radii[6] * uL)
            push!(ArrayR7_temp, radii[7] * uL)
            push!(ArrayR8_temp, radii[8] * uL)
            push!(ArrayR9_temp, radii[9] * uL)

            total_halo_mass = sum(rho) * unit_cell_volumn
            ax_all, ay_all, az_all = grad_central(-Δ..., Φ_all)
            plot_virial && update_virial_terms!(ArrayVirialPotential_temp, ArrayTotalKineticE_temp, ArrayTotalQuantumE_temp, ArrayVirial_temp, ArrayMomentumX_temp, ArrayMomentumY_temp, ArrayMomentumZ_temp, ψ, Φ_all, rho, sqrt_rho, Δ, unit_cell_volumn, uMomentum)

            if average
                if t[i] * time_astro >= average_start_t
                    buffer_ψ2 += rho
                    average_N += 1
                    if average_snapshots
                        save(joinpath(outputdir, "$(title), $(suffix) - timestep_$(@sprintf("%06d", i)).jld2"),
                            "ψ", ψ,
                            "Φ_all", Φ_all,
                        )
                    end
                end
            end

            # Update plot
            if Realtime && iszero(mod(i, StepsBetweenSnapshots))
                if distributed_memory
                    SliceXY[] = dropdims(collect(rho[:, :, rho_max_id[3]]), dims = 3)
                    SliceYZ[] = dropdims(collect(rho[rho_max_id[1], :, :]), dims = 1)
                    SliceXZ[] = dropdims(collect(rho[:, rho_max_id[2], :]), dims = 2)
                else
                    SliceXY[] = rho[:, :, rho_max_id[3]]
                    SliceYZ[] = rho[rho_max_id[1], :, :]
                    SliceXZ[] = rho[:, rho_max_id[2], :]
                end

                r_filter = r_mass_center .<= target_profile_rs * target_fitting_rs_ratio / length_astro

                _profile_r_mean, _profile_ρ_mean, _profile_r_std, _profile_ρ_std = distribution(r_mass_center[r_filter] * uL, collect(rho)[r_filter] * uRho;
                    section = ceil(Int, target_profile_rs * target_fitting_rs_ratio / length_astro / Δ[1]),
                )
                profile_r_mean[] = _profile_r_mean
                profile_ρ_mean[] = _profile_ρ_mean

                append!(ArrayTotalMass[], ArrayTotalMass_temp); empty!(ArrayTotalMass_temp)
                append!(ArrayR[], ArrayR_temp);   empty!(ArrayR_temp)
                append!(ArrayR1[], ArrayR1_temp); empty!(ArrayR1_temp)
                append!(ArrayR2[], ArrayR2_temp); empty!(ArrayR2_temp)
                append!(ArrayR3[], ArrayR3_temp); empty!(ArrayR3_temp)
                append!(ArrayR4[], ArrayR4_temp); empty!(ArrayR4_temp)
                append!(ArrayR5[], ArrayR5_temp); empty!(ArrayR5_temp)
                append!(ArrayR6[], ArrayR6_temp); empty!(ArrayR6_temp)
                append!(ArrayR7[], ArrayR7_temp); empty!(ArrayR7_temp)
                append!(ArrayR8[], ArrayR8_temp); empty!(ArrayR8_temp)
                append!(ArrayR9[], ArrayR9_temp); empty!(ArrayR9_temp)

                ArrayT[] = t[1:i] * uT
                push!(ArrayT_Snap[], t[i] * uT)

                plot_virial && update_virial_visualization!(ArrayVirialPotential, ArrayVirialPotential_temp, ArrayTotalKineticE, ArrayTotalKineticE_temp, ArrayTotalQuantumE, ArrayTotalQuantumE_temp, ArrayVirial, ArrayVirial_temp, ArrayMomentumX, ArrayMomentumX_temp, ArrayMomentumY, ArrayMomentumY_temp, ArrayMomentumZ, ArrayMomentumZ_temp, AxisVirial, ArrayT)
                
                Makie.xlims!(AxisR, 0, ArrayT[][end])
                Makie.ylims!(AxisR, minimum(ArrayR1[]), maximum(ArrayR9[]))

                if dynamic_colorrange
                    ColorRange[] = (0, maximum(abs.(ψ).^2)/10)
                end
                recordframe!(io)
            end

            if extract_dwarf_granule && (t[i] * time_astro >= extract_min_t)
                if extract_mode == :profile
                    current_fit_error = compute_profile_fit_error(r_mass_center, rho, length_astro, Δ, density_astro, profile_config, target_profile_model, uniform_interval)
                elseif extract_mode == :RC
                    current_fit_error = compute_rc_fit_error(r_mass_center, ax_all, ay_all, az_all, xxx, yyy, zzz, rho_max_id, length_astro, Δ, astro_config, df_CO_RC, uniform_interval)
                end
                best_fit_error, best_fit_t, best_fit_beta_star_error, best_fit_beta_star = update_best_fit!(best_fit_error, current_fit_error, t, i, time_astro, best_fit_ψ, ψ, best_fit_ψ_last_t, ψ_last_t, best_fit_Φ_all, Φ_all, rc_config, fig, outputdir, title, suffix, r_mass_center, rho, length_astro)
            end

            update_unicode_progress!(progress, i, t, unicode_plot, distributed_memory, rho, rho_max_id, Realtime, StepsBetweenSnapshots, r_target, ρ_halo_target, _profile_r_mean, _profile_ρ_mean, best_fit_t, best_fit_error, current_fit_error, best_fit_beta_star_error, best_fit_beta_star, current_beta_star, unicode_heatmap_width, Xmax, uT, uL, Nx, Δ)

            if need_to_interrupt(outputdir, remove = true)
                breakflag = true
                break
            end
        end
    end

    unicode_plot && println("\n"^(div(unicode_heatmap_width,2)))

    dfProp = save_property_dataframe(ArrayT, ArrayR, ArrayR1, ArrayR2, ArrayR3, ArrayR4, ArrayR5, ArrayR6, ArrayR7, ArrayR8, ArrayR9, ArrayTotalMass, plot_virial, ArrayVirialPotential, ArrayTotalKineticE, ArrayTotalQuantumE, ArrayVirial, ArrayMomentumX, ArrayMomentumY, ArrayMomentumZ, outputdir, title, suffix)
    averaged_ψ2, averaged_a_all = compute_averaged_fields(average, buffer_ψ2, average_N, baryon_mode, a_all, Φ_b, Δ, Nx, Ny, Nz, gpu, GPU, CPU, Periodic, fft_poisson, grad_central)

    if baryon_mode == :ignored
        # z=0, y=0 plane
        dfAcc = DataFrame(
            :r => vec(r[:, :, div(end,2)]) .* uL,
            :a_all => vec(a_all) .* uAcc,
            :a_all_averaged => vec(averaged_a_all) .* uAcc,
        )
        chi2RC = nothing
    else
        # z=0, y=0 plane
        dfAcc = DataFrame(
            :r => rMOND[:],
            :a_b => a_b[:] .* uAcc,
            :a_all => a_all[:] .* uAcc,
            :a_all_averaged => averaged_a_all[:] .* uAcc,
            :a_mond => a_mond[:] .* uAcc,
        )

        mass_fix_ratio = 1
        figRC, chi2RC, mass_fix_ratio = plot_MW_RC_SPE(dfAcc; best_fit_halo_mass)
        Makie.save(joinpath(outputdir, "$(title), $(suffix) - RC.png"), figRC)

        @info "taking average"
        figRC, _chi2RC, _mass_fix_ratio = plot_MW_RC_SPE(dfAcc; average = true, best_fit_halo_mass)
        Makie.save(joinpath(outputdir, "$(title), $(suffix) - RC averaged.png"), figRC)

        figRAR = compute_RAR(dfAcc; minR, maxR, plotMaxR = massRadius, zoom_a_max, mass_fix_ratio)
        Makie.save(joinpath(outputdir, "$(title), $(suffix) - RAR.png"), figRAR)

        if average
            figRAR = compute_RAR(dfAcc; minR, maxR, plotMaxR = massRadius, zoom_a_max, mass_fix_ratio, average = true)
            Makie.save(joinpath(outputdir, "$(title), $(suffix) - RAR averaged.png"), figRAR)
        end
    end
    CSV.write(joinpath(outputdir, "$(title), $(suffix) - acc.csv"), dfAcc)

    save_phi && save_evolution_results(ψ, Φ_all, ψ_last_t, average, averaged_ψ2, outputdir, title, suffix)

    Makie.save(joinpath(outputdir, "$(title), $(suffix) - Overview Prop.png"), fig)

    if baryon_mode == :ignored
        figMOND, MOND_errorrel = nothing, nothing
    else
        figMOND, MOND_errorrel = plotMOND(ax_all, ay_all, az_all, ax_b, ay_b, az_b, a0, r, length_astro, acc_astro, minR, maxR, outputdir, title, suffix, section; filename = "$(title)_Prop")
    end
    @info "Files saved to folder: $(outputdir)"

    if extract_dwarf_granule
        println()
        @info "best_fit_error = $(best_fit_error)"
        @info "best_fit_t = $(best_fit_t)"

        ## Finally save data
        if isinf(best_fit_error)
            @warn "Dwarf galaxy: Optimization failed!"
        else
            save(joinpath(outputdir, "$(title), $(suffix) - Prop best fit.jld2"), "ψ", best_fit_ψ, "Φ_all", best_fit_Φ_all, "ψ_last_t", best_fit_ψ_last_t)

            write(best_fit_time_file, @sprintf("%.4f\n", ustrip(u"Gyr", best_fit_t)))
            close(best_fit_time_file)
        end
    end

    if Realtime
        display(fig)
    end
    return ψ, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2
end


function test_MW_MOND(;
    model = :MW,
    V = (x,y,z,ψ)->0,
    FDM_mass_ratio = 1.0,
    FDM_radius_ratio = 1.0,
    title = "WaveDM_MOND_MW",
    Xmax = 2.0,
    Ymax = Xmax,
    Zmax = Xmax,
    Tmax = 0.1,
    Nt = 5000,
    Nx = 256,
    Ny = Nx,
    Nz = Nx,
    Np = 5000, # Number of particles of each baryon component
    absorb_coeff = 10.0,
    StepsBetweenSnapshots = 5,
    IC_vel = nothing, # If nothing, generate a new IC; otherwise, use the provided IC
    ρ_baryon = nothing,
    reset_velocity = false, # If true, reset velocity distribution
    Φ_b = nothing,
    ax_b = nothing,
    ay_b = nothing,
    az_b = nothing,
    IC_only = false, # If true, return IC without evolution
    static = false, # If true, no velocity, phases are all zero
    save_IC = true,
    rotational_ratio = 0.0,
    velocity_ratio = 1.0,
    velocity_falling = false, # if false, use use random direction; if true, all velocities point to zero point
    outputdir = joinpath(@__DIR__, "output/MOND"),
    massRadius = 50u"kpc",
    bulk_perturb = true, # If true, pooling the initial conditions
    bulk_size = 8,
    bulk_center_size = 0,
    bulk_shift_size = div(bulk_size,2)-1,
    baryon_β = 0.65,
    baryon_ρ0 = 1.18e6u"Msun/kpc^3",
    baryon_r0 = 100u"kpc",
    halo_β = 0.5,
    halo_ρ0 = 3e5u"Msun/kpc^3",
    halo_r0 = 600u"kpc",
    halo_α = 1.0,
    halo_γ = 1.0,
    halo_Q = 1.0,
    baryon_fraction_limit = 1,
    GravitySolver = Tree(),
    boundary = Periodic(),
    SofteningLength = 1.0u"kpc",
    baryon_mode = :mesh,
    gpu = true,

    ### dwarfs
    stellar_TotalMass = NaN,
    stellar_ScaleRadius = NaN,
    thickness_ratio_stellar = NaN,

    gases_TotalMass = NaN,
    gases_ScaleRadius = NaN,
    thickness_ratio_gases = NaN,

    mₐ = 0.2 * 1e-22 * 1.783e-36u"kg",
    Ωₘ₀ = 0.31,
    aₛ = +8.29e-60u"fm", # FDM

    length_astro = uconvert(u"kpc", (8 * π * (C.h/2/π)^2 / (3 * mₐ^2 * C.H^2 * Ωₘ₀))^0.25),
    time_astro = uconvert(u"Gyr", (3 * C.H^2 * Ωₘ₀ / (8 * π))^-0.5),
    mass_astro = uconvert(u"Msun", (3 * C.H^2 * Ωₘ₀ / (8 * π))^0.25 * (C.h/2/π)^1.5 / (mₐ^1.5 * C.G)),
    density_astro = uconvert(u"Msun/kpc^3", mass_astro/length_astro^3),
    acc_astro = uconvert(u"m/s^2", length_astro / time_astro^2),
    velocity_astro = uconvert(u"m/s", length_astro / time_astro),
    potential_astro = uconvert(u"kpc^2/Gyr^2", C.G * mass_astro / length_astro),

    h_astro = C.h / length_astro^2 / mass_astro * time_astro |> upreferred,
    aₛ_astro = aₛ / length_astro |> upreferred,
    mₐ_astro = mₐ / mass_astro |> upreferred,
    c_astro = C.c / length_astro * time_astro |> upreferred,
    κ_astro = 4π * ħ * aₛ / (time_astro * mₐ^2 * C.G) |> upreferred,

    G0 = C.G * mass_astro * time_astro^2 / length_astro^3,
    a0_astro = C.ACC0 / length_astro * time_astro^2,

    distributed_memory = false,
    pids = workers(),
    kw...
)

    @time if isnothing(IC_vel)
        @info "Initializing grid"
        x, y, z, Δ, unit_cell_volumn = setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)
        dt = Tmax / Nt
        t = collect(LinRange(0, Tmax, Nt))
        DA = distributed_memory ? DArray : collect
        oneMatrix = ones(Nx, Ny, Nz)
        xxx, yyy, zzz, rrr = setup_coordinates(x, y, z, Nx, Ny, Nz, oneMatrix; DA)
        # RRR = sqrt.(xxx.^2 + yyy.^2)
        
        r_in_range = collect(rrr .< upreferred(massRadius / length_astro))

        @info "Sampling IC density"
        
        # Define unit conversion factors (needed for AstroUnitsConfig)
        uT = ustrip(u"Gyr", time_astro)
        uL = ustrip(u"kpc", length_astro)
        uVel = ustrip(u"km/s", velocity_astro)
        uAcc = ustrip(u"m/s^2", acc_astro)
        uRho = ustrip(u"Msun/kpc^3", density_astro)
        uMomentum = ustrip(u"Msun*kpc/Gyr", mass_astro * velocity_astro)
        
        # Create astrophysical units configuration
        astro_config = AstroUnitsConfig(length_astro, time_astro, mass_astro, density_astro, acc_astro, velocity_astro, potential_astro, 
                                        uT, uL, uVel, uAcc, uRho, uMomentum, h_astro, aₛ_astro, mₐ_astro, c_astro, κ_astro, G0, a0_astro)
        
        # Create simulation grid
        grid = SimulationGrid(Xmax, Ymax, Zmax, Nx, Ny, Nz, Δ, x, y, z, xxx, yyy, zzz, rrr, oneMatrix, unit_cell_volumn)
        
        # Create initial conditions configuration
        init_config = InitialConditionsConfig(model, baryon_mode, Np, pids, bulk_perturb, bulk_size, bulk_shift_size, bulk_center_size, 
                                            reset_velocity, static, FDM_mass_ratio, FDM_radius_ratio, GravitySolver, SofteningLength)
        
        # Create density profile configuration - use original physical quantities with units
        density_config = DensityProfileConfig(
            baryon_β, baryon_ρ0, baryon_r0, halo_β, halo_ρ0, halo_r0, halo_α, halo_γ, halo_Q, stellar_TotalMass, stellar_ScaleRadius, thickness_ratio_stellar, gases_TotalMass, gases_ScaleRadius, thickness_ratio_gases
        )
        
        # Call generate_initial_conditions with new struct parameters
        ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon = generate_initial_conditions(init_config, grid, density_config, astro_config)
        
        total_mass_halo_IC = sum(ρ_halo[r_in_range]) * prod(Δ) * mass_astro
        
        @info "Total mass of halo: $(total_mass_halo_IC)"
        if baryon_mode != :ignored
            baryon_ratio_IC = total_mass_baryon / (total_mass_baryon + total_mass_halo_IC)
            @info "Total mass of baryon: $(total_mass_baryon)"
            @info "Baryon ratio: $(@sprintf("%.2f", baryon_ratio_IC * 100)) %"

            if baryon_ratio_IC > baryon_fraction_limit
                return nothing
            end
        else
            baryon_ratio_IC = 0
        end

        @info "Computing acc"
        @time Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM = compute_acceleration_field(ρ_halo, Δ, boundary, gpu, xxx, yyy, zzz, unit_cell_volumn, mass_astro, SofteningLength, potential_astro, pids, length_astro)
        
        @info "Setting velocities"
        v, Φ_all, ax_all, ay_all, az_all = generate_velocity_field(Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM, Φ_b, ax_b, ay_b, az_b, baryon_mode, xxx, yyy, zzz, Δ, velocity_falling, rotational_ratio, velocity_ratio)
        vx, vy, vz = adjust_velocity_field(v, ρ_halo, bulk_perturb, bulk_size, bulk_shift_size, bulk_center_size)
        
        @info "Solving vector equation of velocity field"
        θ = solve_vector_equation(collect(vx), collect(vy), collect(vz), Δ...)
        phase = exp.(im*θ)
        IC_vel = sqrt.(ρ_halo) .* phase

        r_in_range = Φ_WaveDM = ax_WaveDM = ay_WaveDM = az_WaveDM = v = Φ_all = ax_all = ay_all = az_all = vx = vy = vz = θ = phase = nothing # release memory
        oneMatrix = xxx = yyy = zzz = rrr = nothing
        if IC_only
            return IC_vel
        end
    elseif IC_vel isa String
        if baryon_mode == :ignored
            IC_vel = load(IC_vel,
                "ψ",
            );
        else
            IC_vel, Φ_b, ax_b, ay_b, az_b = load(IC_vel,
                "ψ", "Φ_b", "ax_b", "ay_b", "az_b",
            );
        end
        save_IC = false
        ρ_baryon = nothing
        total_mass_baryon = 0.0u"Msun"
        baryon_ratio_IC = 0

        if reset_velocity
            x, y, z, Δ, unit_cell_volumn = setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)
            t = collect(LinRange(0, Tmax, Nt))
            dt = Tmax / Nt
            DA = distributed_memory ? DArray : collect
            oneMatrix = ones(Nx, Ny, Nz)
            xxx, yyy, zzz, rrr = setup_coordinates(x, y, z, Nx, Ny, Nz, oneMatrix; DA)
            
            ρ_halo = abs.(IC_vel).^2
            Φ_WaveDM = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], ρ_halo, Periodic(), gpu ? GPU() : CPU()))

            ax_WaveDM, ay_WaveDM, az_WaveDM = grad_central(-Δ..., Φ_WaveDM)

            @info "Setting velocities"
            v, Φ_all, ax_all, ay_all, az_all = generate_velocity_field(Φ_WaveDM, ax_WaveDM, ay_WaveDM, az_WaveDM, Φ_b, ax_b, ay_b, az_b, baryon_mode, xxx, yyy, zzz, Δ, velocity_falling, rotational_ratio, velocity_ratio)
            vx, vy, vz = adjust_velocity_field(v, ρ_halo, bulk_perturb, bulk_size, bulk_shift_size, bulk_center_size)
                
            @info "Solving vector equation of velocity field"
            θ = solve_vector_equation(vx, vy, vz, Δ...)
            phase = exp.(im*θ)
            IC_vel = sqrt.(ρ_halo) .* phase

            Φ_WaveDM = ax_WaveDM = ay_WaveDM = az_WaveDM = v = Φ_all = ax_all = ay_all = az_all = vx = vy = vz = θ = phase = nothing # release memory
            oneMatrix = xxx = yyy = zzz = rrr = nothing
        end
    else  # Directly use the IC_vel variable
        ρ_baryon = nothing
        total_mass_baryon = 0.0u"Msun"
    end


    @info "Start SPE simulation"
    return SPE3D_MOND(;
        Xmax, Ymax, Zmax, Tmax, Nx, Ny, Nz, Nt,
        limsX = (-0.3, 0.3),
        limsY = (-0.3, 0.3),
        limsZ = (-0.3, 0.3),
        StepsBetweenSnapshots, absorb_coeff, IC = IC_vel, save_IC, baryon = ρ_baryon, baryon_mode, baryon_potential = Φ_b,
        ax_b, ay_b, az_b, boundary, V, title, outputdir, plotWaveDM = true,
        FDM_mass_ratio, FDM_radius_ratio, rotational_ratio, massRadius, length_astro, time_astro, mass_astro, density_astro, acc_astro, velocity_astro, potential_astro,
        h_astro, aₛ_astro, mₐ_astro, c_astro, κ_astro, G0, a0_astro, gpu, GravitySolver, distributed_memory, pids, kw...
    )
end
