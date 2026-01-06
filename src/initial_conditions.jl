# MOND initial conditions module

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
