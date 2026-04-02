# SPE tidal fields module

"""
$(TYPEDSIGNATURES)

Setup Milky Way tidal field interpolation.
This function encapsulates the MW tidal field setup code from SPE3D_waveDM.
"""
function setup_mw_tidal_field(MW_pot, MW_pot_Xmax, MW_pot_Ymax, MW_pot_Zmax, MW_pot_N,
    length_astro, potential_astro, spl_pot, sim_force_baryon, SofteningLength)
    if isnothing(MW_pot)
        @info "MW tidal force: computing potential on interpolation grid"
        MW_x = Vector{Float64}(LinRange(-MW_pot_Xmax, MW_pot_Xmax, MW_pot_N)) ./ length_astro
        MW_y = Vector{Float64}(LinRange(-MW_pot_Ymax, MW_pot_Ymax, MW_pot_N)) ./ length_astro
        MW_z = Vector{Float64}(LinRange(-MW_pot_Zmax, MW_pot_Zmax, MW_pot_N)) ./ length_astro

        MW_pos = PVector.(Base.Iterators.product(MW_x*length_astro, MW_y*length_astro, MW_z*length_astro))

        @time MW_Phi_DM = spl_pot.(ustrip.(u"kpc", norm.(MW_pos))) * u"kpc^2/Gyr^2" / potential_astro
        if isnothing(sim_force_baryon)
            MW_Phi = MW_Phi_DM
        else
            @info "MW tidal force: computing potentials from baryons"
            @time MW_Phi_baryon = compute_potential(sim_force_baryon, MW_pos, SofteningLength, sim_force_baryon.config.solver.grav, CPU()) / potential_astro
            MW_Phi = MW_Phi_DM + MW_Phi_baryon
        end
    elseif MW_pot isa AbstractString
        @info "Loading MW tidal potential field from file: $(MW_pot)"
        MW_x, MW_y, MW_z = load(MW_pot, "MW_x", "MW_y", "MW_z") ./ length_astro
        MW_Phi = load(MW_pot, "MW_Phi") / potential_astro
    end

    MW_grid = RectangleGrid(MW_x, MW_y, MW_z)

    return MW_grid, MW_Phi, MW_x, MW_y, MW_z
end

"""
$(TYPEDSIGNATURES)

Compute tidal potential at current simulation time using config objects.
Accesses config fields directly without internal unpacking.
"""
function add_tidal_potential!(
    Φ_all,
    config_tidal::TidalFieldConfig,
    grid::SimulationGrid,
    config_gravity::GravityConfig,
    i::Int,
    t::Vector{<:Real}
)
    id_t = findfirstvalue(ustrip(u"Gyr", config_tidal.tidal_lookback_time) .- config_tidal.df_traj.time, t[i] * config_tidal.uT)
    if isnothing(id_t)
        @warn "Current time $(t[i]*time_astro): setting lookback time = 0 Gyr"
        id_t = 1
    end

    traj_pos = PVector{Float64}[]
    if config_tidal.LMC_tidal_field || !config_tidal.MW_tidal_interpolate
        traj_center = PVector(config_tidal.df_traj.x[id_t], config_tidal.df_traj.y[id_t], config_tidal.df_traj.z[id_t]) * u"kpc"
        traj_pos = PVector.(grid.xxx * config_tidal.length_astro, grid.yyy * config_tidal.length_astro, grid.zzz * config_tidal.length_astro) .+ traj_center
    end

    if config_tidal.MW_tidal_interpolate
        pot_tidal = GridInterpolations.interpolate.(Ref(config_tidal.MW_grid), Ref(config_tidal.MW_Phi), collect.(zip(grid.xxx .+ config_tidal.df_traj.x[id_t]/config_tidal.uL, grid.yyy .+ config_tidal.df_traj.y[id_t]/config_tidal.uL, grid.zzz .+ config_tidal.df_traj.z[id_t]/config_tidal.uL)))
    else
        traj_r = norm.(traj_pos)
        pot_DM = config_tidal.spl_pot.(ustrip.(u"kpc", traj_r)) * u"kpc^2/Gyr^2"
        if isnothing(config_tidal.sim_force_baryon)
            pot_tidal = pot_DM / config_gravity.mass_astro
        else
            pot_b = compute_potential(config_tidal.sim_force_baryon, traj_pos, config_gravity.SofteningLength, config_tidal.sim_force_baryon.config.solver.grav, CPU())
            pot_tidal = (pot_b + pot_DM) / config_gravity.mass_astro
        end
    end

    if config_tidal.LMC_tidal_field
        id_t_LMC = findfirstvalue(ustrip(u"Gyr", config_tidal.tidal_lookback_time) .- config_tidal.df_traj_LMC.time, t[i] * config_tidal.uT)
        if ustrip(u"Gyr", config_tidal.tidal_lookback_time) <= t[i] * config_tidal.uT
            id_t_LMC = 1
        end
        if config_gravity.GravitySolver isa DirectSum
            config_tidal.sim_traj_LMC.simdata.Pos .= config_tidal.particles_LMC.Pos .+ PVector(config_tidal.df_traj_LMC.x[id_t_LMC], config_tidal.df_traj_LMC.y[id_t_LMC], config_tidal.df_traj_LMC.z[id_t_LMC]) * u"kpc"
        elseif config_gravity.GravitySolver isa Tree
            config_tidal.sim_traj_LMC.simdata.tree.data.Pos .= config_tidal.particles_LMC.Pos .+ PVector(config_tidal.df_traj_LMC.x[id_t_LMC], config_tidal.df_traj_LMC.y[id_t_LMC], config_tidal.df_traj_LMC.z[id_t_LMC]) * u"kpc"
            AstroNbodySim.rebuild_tree(config_tidal.sim_traj_LMC)
        end
        pot_LMC = AstroNbodySim.compute_potential(config_tidal.sim_traj_LMC, traj_pos, config_gravity.SofteningLength, config_gravity.GravitySolver, CPU()) / config_gravity.mass_astro
        pot_tidal = pot_tidal + pot_LMC
    end

    pot_tidal .-= pot_tidal[grid.rho_max_id]
    cancel_field_gradient_at_center!(pot_tidal, grid.rho_max_id, grid.oneMatrix, grid.Nx, grid.Ny, grid.Nz)

    Φ_all += pot_tidal
end

"""
$(TYPEDSIGNATURES)

Cancel out gradient of a field at a specified center point.
This is a generic function that can be used for both potential fields and tidal fields.
"""
function cancel_field_gradient_at_center!(field, center_id, oneMatrix, Nx, Ny, Nz)
    dp_dx = (field[center_id[1]+1, center_id[2],   center_id[3]]   - field[center_id[1]-1, center_id[2],   center_id[3]])/2
    dp_dy = (field[center_id[1],   center_id[2]+1, center_id[3]]   - field[center_id[1],   center_id[2]-1, center_id[3]])/2
    dp_dz = (field[center_id[1],   center_id[2],   center_id[3]+1] - field[center_id[1],   center_id[2],   center_id[3]-1])/2

    field -= oneMatrix .* (collect(1:Nx) .- div(Nx,2)) * dp_dx
    field -= oneMatrix .* (collect(1:Ny) .- div(Ny,2))' * dp_dy
    field -= oneMatrix .* reshape(collect(1:Nz) .- div(Nz,2), 1, 1, Nz) * dp_dz

    return field
end