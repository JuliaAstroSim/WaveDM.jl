# Milky Way IC generation functions

"""
$(TYPEDSIGNATURES)

Compute total Milky Way density at position (x, y, z).
"""
function density_all_MW(x, y, z)
    ρ_bulge, ρ_discs, ρ_gas, ρ_halo = GalacticDynamics.milkyway_zhu2023how()
    R = sqrt(x^2 + y^2)
    r = sqrt(x^2 + y^2 + z^2)
    return ρ_bulge(R,z) + ρ_discs(R,z) + ρ_gas(R,z) + ρ_halo(r)
end

"""
$(TYPEDSIGNATURES)

Compute baryonic Milky Way density at position (x, y, z).
"""
function density_baryon_MW(x, y, z)
    ρ_bulge, ρ_discs, ρ_gas, ρ_halo = GalacticDynamics.milkyway_zhu2023how()
    R = sqrt(x^2 + y^2)
    return ρ_bulge(R,z) + ρ_discs(R,z) + ρ_gas(R,z)
end

"""
$(TYPEDSIGNATURES)

Generate Milky Way initial conditions.
"""
function generate_milkyway_initial_conditions(grid::SimulationGrid, config_IC::InitialConditionsConfig, config_units::AstroUnitsConfig, config_device::DeviceConfig, boundary)
    length_astro = config_units.length_astro
    density_astro = config_units.density_astro
    potential_astro = config_units.potential_astro
    acc_astro = config_units.acc_astro

    model_halo = Zhao(1.55e7u"Msun/kpc^3" * config_IC.FDM_mass_ratio, 11.75u"kpc" * config_IC.FDM_radius_ratio, 1.19, 2.95, 0.95)
    ρ_halo = sampling_density.(grid.r, model_halo, length_astro, density_astro) |> collect
    baryon_particles = nothing
    if config_IC.baryon_mode == :mesh
        ρ_baryon = density_baryon_MW.(grid.xxx*length_astro, grid.yyy*length_astro, grid.zzz*length_astro)/density_astro
        Φ_b = 4π * parallel_poisson(grid.Δ, [grid.Nx-1, grid.Ny-1, grid.Nz-1], ρ_baryon, boundary, config_device)
        ax_b, ay_b, az_b = grad_central(-grid.Δ..., Φ_b)
        total_mass_baryon = sum(ρ_baryon) * grid.unit_cell_volumn * density_astro
    elseif config_IC.baryon_mode == :particles_static || config_IC.baryon_mode == :particles_dynamic
        pos = PVector.(grid.xxx * length_astro, grid.yyy * length_astro, grid.zzz * length_astro)
        baryon_particles = generate_milkyway_baryon_particles(config_IC.Np, config_IC)

        sim_force_baryon = Simulation(baryon_particles; GravitySolver = config_IC.GravitySolver, pids = config_IC.pids)

        @info "Computing baryonic potentials and forces with $(traitstring(config_IC.GravitySolver)) solver"
        @time Φ_b = compute_potential(sim_force_baryon, pos, config_IC.SofteningLength, config_IC.GravitySolver, CPU()) ./ potential_astro
        @time acc_b = StructArray(compute_force(sim_force_baryon, pos, config_IC.SofteningLength, config_IC.GravitySolver, CPU()))
        ax_b = upreferred.(acc_b.x ./ acc_astro)
        ay_b = upreferred.(acc_b.y ./ acc_astro)
        az_b = upreferred.(acc_b.z ./ acc_astro)

        total_mass_baryon = sum(baryon_particles.Mass)
        ρ_baryon = nothing
    elseif config_IC.baryon_mode == :ignored
        ρ_baryon = Φ_b = ax_b = ay_b = az_b = baryon_particles = nothing
        total_mass_baryon = 0.0u"Msun"
    end

    return ρ_halo, ρ_baryon, Φ_b, ax_b, ay_b, az_b, total_mass_baryon, baryon_particles
end

# Generate Milky Way baryon particles
function generate_milkyway_baryon_particles(Np::Int, config_IC)
    TotalMass_bulge = 8.5708e9u"Msun"
    TotalMass_thin = uconvert(u"Msun", 2pi * 1003.12u"Msun/pc^2" * (2.42u"kpc")^2) # 3.691165259383738e10 M⊙
    TotalMass_thick = uconvert(u"Msun", 2pi * 167.93u"Msun/pc^2" * (3.17u"kpc")^2) # 1.0602949202938915e10 M⊙
    TotalMass_HI = 1.0674e10u"Msun"
    TotalMass_HII = 1.2303e9u"Msun"
    TotalMass_baryons = TotalMass_bulge + TotalMass_thin + TotalMass_thick + TotalMass_HI + TotalMass_HII

    NumSamples_bulge = ceil(Int, TotalMass_bulge/TotalMass_baryons * Np)
    NumSamples_thin  = ceil(Int, TotalMass_thin/TotalMass_baryons * Np)
    NumSamples_thick = ceil(Int, TotalMass_thick/TotalMass_baryons * Np)
    NumSamples_HI    = ceil(Int, TotalMass_HI/TotalMass_baryons * Np)
    NumSamples_HII = Np - NumSamples_bulge - NumSamples_thin - NumSamples_thick - NumSamples_HI

    @info "NumSamples of bulge: $(NumSamples_bulge)"
    @info "NumSamples of thin:  $(NumSamples_thin)"
    @info "NumSamples of thick: $(NumSamples_thick)"
    @info "NumSamples of HI:    $(NumSamples_HI)"
    @info "NumSamples of HII:   $(NumSamples_HII)"

    if config_IC.MW_disk_RC
        @info "Setting rotation velocities of stellar and gaseous disks (Eilers et al. 2019)"
        RC_MW = load_MW_RC_Eilers2019()
        RotationCurve = ([0.0; RC_MW.r] * u"kpc", [0.0; RC_MW.vel] * u"km/s"); # add a zero point
    else
        RotationCurve = nothing
    end

    particles_bulge = generate(AstroIC.Bulge(;
        collection = BULGE,
        NumSamples = NumSamples_bulge,
        TotalMass = TotalMass_bulge,
        ScaleRadius = 0.075u"kpc",
        CutRadius = 2.1u"kpc",
        q = 0.5,
        α = 1.8,
    ))

    particles_stellar_thin = generate(AstroIC.ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_thin,
        TotalMass = TotalMass_thin,
        ScaleRadius = 2.42u"kpc",
        ScaleHeight = 0.3u"kpc",
    ); RotationCurve, rotational_ratio = config_IC.MW_disk_RC)

    particles_stellar_thick = generate(AstroIC.ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_thick,
        TotalMass = TotalMass_thick,
        ScaleRadius = 3.17u"kpc",
        ScaleHeight = 0.9u"kpc",
    ); RotationCurve, rotational_ratio = config_IC.MW_disk_RC)

    particles_gas_HI = generate(AstroIC.ExponentialDisc(;
        collection = GAS,
        NumSamples = NumSamples_HI,
        TotalMass = TotalMass_HI,
        ScaleRadius = 7.0u"kpc",
        ScaleHeight = 0.085u"kpc",
        HoleRadius = 4.0u"kpc",
    ); RotationCurve, rotational_ratio = config_IC.MW_disk_RC)

    particles_gas_HII = generate(AstroIC.ExponentialDisc(;
        collection = GAS,
        NumSamples = NumSamples_HII,
        TotalMass = TotalMass_HII,
        ScaleRadius = 1.5u"kpc",
        ScaleHeight = 0.045u"kpc",
        HoleRadius = 12.0u"kpc",
    ); RotationCurve, rotational_ratio = config_IC.MW_disk_RC)

    particles = [particles_bulge; particles_stellar_thin; particles_stellar_thick; particles_gas_HI; particles_gas_HII]

    return particles
end