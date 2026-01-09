# SPE visualization module

"""
$(TYPEDSIGNATURES)

Setup initial visualization figure with multiple panels.
"""
function setup_visualization(grid::SimulationGrid, t::Vector, vis_config::VisualizationConfig, data_config::VisualizationData, astro_config::AstroUnitsConfig, distributed_memory::Bool)::Tuple{Figure, Observable{Vector{Any}}, Observable{Vector{Any}}, Axis, Axis, Axis, Observable{Matrix{Float64}}, Observable{Matrix{Float64}}, Observable{Matrix{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Vector{Float64}}, Observable{Tuple{Float64, Float64}}}
    Xmax = grid.Xmax
    Δ = grid.Δ
    x = grid.x
    y = grid.y
    z = grid.z
    r = grid.r
    rho = data_config.rho
    rho_max_id = data_config.rho_max_id
    total_halo_mass = data_config.total_halo_mass
    radii = data_config.radii
    uT = astro_config.uT
    uL = astro_config.uL
    title = vis_config.title
    suffix = vis_config.suffix
    size = vis_config.size
    fig = Figure(; size)
    
    ArrayT = Observable([t[1] * uT])
    ArrayT_Snap = Observable([t[1] * uT])

    AxisXY = Makie.Axis(fig[1,1];
        title = "", xlabel = "x [kpc]", ylabel = "y [kpc]", aspect = 1,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    AxisYZ = Makie.Axis(fig[1,2];
        title = "", xlabel = "y [kpc]", ylabel = "z [kpc]", aspect = 1,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    AxisXZ = Makie.Axis(fig[1,3];
        title = "", xlabel = "x [kpc]", ylabel = "z [kpc]", aspect = 1,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    AxisR = Makie.Axis(fig[2,1];
        xlabel = "t [Gyr]",
        ylabel = "Lagrange radii [kpc]",
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    
    AxisVirial = Makie.Axis(fig[2,2];
        xlabel = "t [Gyr]", ylabel = "virial",
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    AxisDensityProfile = Makie.Axis(fig[2,3];
        xlabel = "log10(r [kpc])", ylabel = "log10(ρ [Msun/kpc^3])", xscale = log10, yscale = log10,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    Makie.xlims!(AxisDensityProfile, 0.5*Δ[1]*uL, Xmax*uL)
    Makie.ylims!(AxisDensityProfile, 1e3, 1e11)
    
    if distributed_memory
        SliceXY = Observable(dropdims(collect(rho[:, :, rho_max_id[3]]), dims = 3))
        SliceYZ = Observable(dropdims(collect(rho[rho_max_id[1], :, :]), dims = 1))
        SliceXZ = Observable(dropdims(collect(rho[:, rho_max_id[2], :]), dims = 2))
    else
        SliceXY = Observable(rho[:, :, rho_max_id[3]])
        SliceYZ = Observable(rho[rho_max_id[1], :, :])
        SliceXZ = Observable(rho[:, rho_max_id[2], :])
    end
    
    ArrayTotalMass = Observable([isinf(total_halo_mass) ? 0 : total_halo_mass])
    # Fix method ambiguity by explicitly calculating weighted std
    r_slice = collect(r[:,:,div(end,2)])[:]
    rho_slice = collect(rho[:,:,div(end,2)])[:]
    w = aweights(rho_slice)
    # Calculate mean manually to avoid sum ambiguity
    mean_val = sum(r_slice .* w) / sum(w)
    var_val = sum(w .* (r_slice .- mean_val).^2) / sum(w)
    std_val = sqrt(var_val)
    ArrayR = Observable([std_val * uL]) #TODO this is strongly affected by rho
    ArrayR1 = Observable([radii[1] * uL])
    ArrayR2 = Observable([radii[2] * uL])
    ArrayR3 = Observable([radii[3] * uL])
    ArrayR4 = Observable([radii[4] * uL])
    ArrayR5 = Observable([radii[5] * uL])
    ArrayR6 = Observable([radii[6] * uL])
    ArrayR7 = Observable([radii[7] * uL])
    ArrayR8 = Observable([radii[8] * uL])
    ArrayR9 = Observable([radii[9] * uL])

    ColorRange = Observable((0, maximum(rho)/10))
    Makie.heatmap!(AxisXY, x * uL, y * uL, SliceXY; colorrange = ColorRange)
    Makie.heatmap!(AxisYZ, y * uL, z * uL, SliceYZ; colorrange = ColorRange)
    hmXZ = Makie.heatmap!(AxisXZ, x * uL, z * uL, SliceXZ; colorrange = ColorRange)
    Makie.lines!(AxisR, ArrayT, ArrayR, label= "std")
    Makie.lines!(AxisR, ArrayT, ArrayR1, label= "Lr1")
    Makie.lines!(AxisR, ArrayT, ArrayR2, label= "Lr2")
    Makie.lines!(AxisR, ArrayT, ArrayR3, label= "Lr3")
    Makie.lines!(AxisR, ArrayT, ArrayR4, label= "Lr4")
    Makie.lines!(AxisR, ArrayT, ArrayR5, label= "Lr5")
    Makie.lines!(AxisR, ArrayT, ArrayR6, label= "Lr6")
    Makie.lines!(AxisR, ArrayT, ArrayR7, label= "Lr7")
    Makie.lines!(AxisR, ArrayT, ArrayR8, label= "Lr8")
    Makie.lines!(AxisR, ArrayT, ArrayR9, label= "Lr9")
    axislegend(AxisR)

    supertitle = Label(fig[0, :], "$(title), $(suffix)")
    return (
        fig,
        ArrayT,
        ArrayT_Snap,
        AxisR,
        AxisVirial,
        AxisDensityProfile,
        SliceXY,
        SliceYZ,
        SliceXZ,
        ArrayTotalMass,
        ArrayR,
        ArrayR1,
        ArrayR2,
        ArrayR3,
        ArrayR4,
        ArrayR5,
        ArrayR6,
        ArrayR7,
        ArrayR8,
        ArrayR9,
        ColorRange,
    )
end

#TODO
"""
$(TYPEDSIGNATURES)

Setup visualization for virial theorem quantities.
"""
function setup_virial_visualization(ψ, Φ_all, rho, sqrt_rho, Δ, unit_cell_volumn, mass_astro, velocity_astro, ArrayT, AxisVirial)
    # Compute initial virial quantities
    ψx, ψy, ψz = grad_central(Δ..., collect(ψ))
    ψc = conj.(ψ)
    ψcx, ψcy, ψcz = grad_central(Δ..., collect(ψc))
    ρvx = real.((ψc .* ψx - ψ .* ψcx) / (2im))  #! slow TODO
    ρvy = real.((ψc .* ψy - ψ .* ψcy) / (2im))
    ρvz = real.((ψc .* ψz - ψ .* ψcz) / (2im))
    # ψx = ψy = ψz = ψc = ψcx = ψcy = ψcz = nothing # release memory

    # Compute momenta
    MomentumX = sum(ρvx) * unit_cell_volumn * mass_astro * velocity_astro
    MomentumY = sum(ρvy) * unit_cell_volumn * mass_astro * velocity_astro
    MomentumZ = sum(ρvz) * unit_cell_volumn * mass_astro * velocity_astro
    ArrayMomentumX = Observable([MomentumX])
    ArrayMomentumY = Observable([MomentumY])
    ArrayMomentumZ = Observable([MomentumZ])
    
    # Compute energies and momenta
    VirialPotential = -0.5*sum(rho[2:end-1,2:end-1,2:end-1] .* Φ_all[2:end-1,2:end-1,2:end-1]) * unit_cell_volumn
    TotalKineticE = 0.5 * sum((ρvx[2:end-1,2:end-1,2:end-1].^2 + ρvy[2:end-1,2:end-1,2:end-1].^2 + ρvz[2:end-1,2:end-1,2:end-1].^2) ./ rho[2:end-1,2:end-1,2:end-1]) * unit_cell_volumn
    sqrt_rho_x, sqrt_rho_y, sqrt_rho_z = grad_central(Δ..., collect(sqrt_rho))
    TotalQuantumE = 0.5 * sum(sqrt_rho_x[2:end-1,2:end-1,2:end-1].^2 + sqrt_rho_y[2:end-1,2:end-1,2:end-1].^2 + sqrt_rho_z[2:end-1,2:end-1,2:end-1].^2) * unit_cell_volumn
    # sqrt_rho_x = sqrt_rho_y = sqrt_rho_z = nothing # release memory
    
    # Setup observables
    ArrayVirialPotential = Observable([VirialPotential])
    ArrayTotalKineticE = Observable([TotalKineticE])
    ArrayTotalQuantumE = Observable([TotalQuantumE])
    ArrayVirial = Observable([2 * TotalKineticE + 2 * TotalQuantumE - VirialPotential]) # Total virial

    Makie.lines!(AxisVirial, ArrayT, ArrayVirialPotential, label = "|V|")
    Makie.lines!(AxisVirial, ArrayT, ArrayTotalKineticE, label = "K")
    Makie.lines!(AxisVirial, ArrayT, ArrayTotalQuantumE, label = "Q")
    Makie.lines!(AxisVirial, ArrayT, ArrayVirial, label = "d²I/dt² = V + 2K + 2Q")
    axislegend(AxisVirial)
    
    return (
        ArrayVirialPotential = ArrayVirialPotential,
        ArrayTotalKineticE = ArrayTotalKineticE,
        ArrayTotalQuantumE = ArrayTotalQuantumE,
        ArrayVirial = ArrayVirial,
        ArrayMomentumX = ArrayMomentumX,
        ArrayMomentumY = ArrayMomentumY,
        ArrayMomentumZ = ArrayMomentumZ,
    )
end

"""
$(TYPEDSIGNATURES)

Update virial theorem quantities during evolution.
"""
function update_virial_terms!(ArrayVirialPotential_temp, ArrayTotalKineticE_temp, ArrayTotalQuantumE_temp, ArrayVirial_temp,
    ArrayMomentumX_temp, ArrayMomentumY_temp, ArrayMomentumZ_temp,
    ψ, Φ_all, rho, sqrt_rho, Δ, unit_cell_volumn, uMomentum)

    # Compute virial quantities
    ψx, ψy, ψz = grad_central(Δ..., collect(ψ))
    ψc = conj.(ψ)
    ψcx, ψcy, ψcz = grad_central(Δ..., collect(ψc))
    ρvx = real.((ψc .* ψx - ψ .* ψcx) / (2im))  #TODO parallel
    ρvy = real.((ψc .* ψy - ψ .* ψcy) / (2im))
    ρvz = real.((ψc .* ψz - ψ .* ψcz) / (2im))

    VirialPotential = -0.5*sum(rho[2:end-1,2:end-1,2:end-1] .* Φ_all[2:end-1,2:end-1,2:end-1]) * unit_cell_volumn
    TotalKineticE = 0.5 * sum((ρvx[2:end-1,2:end-1,2:end-1].^2 + ρvy[2:end-1,2:end-1,2:end-1].^2 + ρvz[2:end-1,2:end-1,2:end-1].^2) ./ rho[2:end-1,2:end-1,2:end-1]) * unit_cell_volumn
    sqrt_rho_x, sqrt_rho_y, sqrt_rho_z = grad_central(Δ..., sqrt_rho) # velocities #TODO: optimize performance
    TotalQuantumE = 0.5 * sum(sqrt_rho_x[2:end-1,2:end-1,2:end-1].^2 + sqrt_rho_y[2:end-1,2:end-1,2:end-1].^2 + sqrt_rho_z[2:end-1,2:end-1,2:end-1].^2) * unit_cell_volumn #TODO: optimize performance
    push!(ArrayVirialPotential_temp, VirialPotential)
    push!(ArrayTotalKineticE_temp, TotalKineticE)
    push!(ArrayTotalQuantumE_temp, TotalQuantumE)
    push!(ArrayVirial_temp, 2 * TotalKineticE + 2 * TotalQuantumE - VirialPotential)

    MomentumX = sum(ρvx) * unit_cell_volumn * uMomentum
    MomentumY = sum(ρvy) * unit_cell_volumn * uMomentum
    MomentumZ = sum(ρvz) * unit_cell_volumn * uMomentum
    push!(ArrayMomentumX_temp, MomentumX)
    push!(ArrayMomentumY_temp, MomentumY)
    push!(ArrayMomentumZ_temp, MomentumZ)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Update visualization for virial theorem quantities during evolution.
"""
function update_virial_visualization!(ArrayVirialPotential, ArrayVirialPotential_temp, ArrayTotalKineticE, ArrayTotalKineticE_temp, ArrayTotalQuantumE, ArrayTotalQuantumE_temp, ArrayVirial, ArrayVirial_temp,
    ArrayMomentumX, ArrayMomentumX_temp, ArrayMomentumY, ArrayMomentumY_temp, ArrayMomentumZ, ArrayMomentumZ_temp, AxisVirial, ArrayT)
    append!(ArrayVirialPotential[], ArrayVirialPotential_temp)
    empty!(ArrayVirialPotential_temp)
    append!(ArrayTotalKineticE[], ArrayTotalKineticE_temp)
    empty!(ArrayTotalKineticE_temp)
    append!(ArrayTotalQuantumE[], ArrayTotalQuantumE_temp)
    empty!(ArrayTotalQuantumE_temp)
    append!(ArrayVirial[], ArrayVirial_temp)
    empty!(ArrayVirial_temp)

    append!(ArrayMomentumX[], ArrayMomentumX_temp)
    empty!(ArrayMomentumX_temp)
    append!(ArrayMomentumY[], ArrayMomentumY_temp)
    empty!(ArrayMomentumY_temp)
    append!(ArrayMomentumZ[], ArrayMomentumZ_temp)
    empty!(ArrayMomentumZ_temp)

    Makie.xlims!(AxisVirial, 0, ArrayT[][end])
    Makie.ylims!(AxisVirial,
        min(minimum(ArrayVirial[]), minimum(ArrayTotalQuantumE[]), minimum(ArrayTotalKineticE[]), minimum(ArrayVirialPotential[])),
        max(maximum(ArrayVirial[]), maximum(ArrayTotalQuantumE[]), maximum(ArrayTotalKineticE[]), maximum(ArrayVirialPotential[])),
    )
end

"""
$(TYPEDSIGNATURES)

Setup visualization for density profile fitting.
"""
function setup_density_profile_visualization!(fig, AxisDensityProfile,
    Δ, rho, r_filter, r_mass_center, target_profile_model, target_profile_error,
    target_profile_ρ0, target_profile_ρ0_u, target_profile_ρ0_d,
    target_profile_rs, target_profile_rs_u, target_profile_rs_d,
    target_profile_α, target_profile_α_u, target_profile_α_d,
    target_profile_β, target_profile_β_u, target_profile_β_d,
    target_profile_γ, target_profile_γ_u, target_profile_γ_d,
    target_fitting_rs_ratio, length_astro, uL, uRho)
    _profile_r_mean, _profile_ρ_mean, _profile_r_std, _profile_ρ_std = distribution(r_mass_center[r_filter] * uL, collect(rho)[r_filter] * uRho;
        # section = ceil(Int, Nx/2*sqrt(3)),
        section = ceil(Int, target_profile_rs * target_fitting_rs_ratio / length_astro / Δ[1]),
    )
    
    profile_r_mean = Observable(deepcopy(_profile_r_mean))
    profile_ρ_mean = Observable(deepcopy(_profile_ρ_mean))

    plots = []
    labels = []

    si1 = Makie.scatter!(AxisDensityProfile, profile_r_mean, profile_ρ_mean)
    li1 = Makie.lines!(AxisDensityProfile, profile_r_mean, profile_ρ_mean)
    push!(plots, [si1, li1])
    push!(labels, "sim")

    r_target = collect(0.05:0.1:ustrip(u"kpc", target_profile_rs * 10))
    if target_profile_model == :dwarf_gNFW
        model_halo_target = gNFW(target_profile_β, target_profile_ρ0, target_profile_rs)
        ρ_halo_target = ustrip.(u"Msun/kpc^3", GalacticDynamics.density.(model_halo_target, r_target*u"kpc"))
    elseif target_profile_model == :dwarf_NFW
        model_halo_target = NFW(target_profile_ρ0, target_profile_rs)
        ρ_halo_target = ustrip.(u"Msun/kpc^3", GalacticDynamics.density.(model_halo_target, r_target*u"kpc"))
    elseif target_profile_model == :dwarf_Zhao
        model_halo_target = Zhao(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ)
        ρ_halo_target = ustrip.(u"Msun/kpc^3", GalacticDynamics.density.(model_halo_target, r_target*u"kpc"))
    elseif target_profile_model == :dwarf_ZhaoQ #TODO consider Q
        model_halo_target = Zhao(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ)
        # model_halo_target = ZhaoQ(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ, target_profile_Q)
        ρ_halo_target = ustrip.(u"Msun/kpc^3", GalacticDynamics.density.(model_halo_target, r_target*u"kpc"))
    end
    lt = Makie.lines!(AxisDensityProfile, r_target, ρ_halo_target, color = :black)
    push!(plots, lt)
    push!(labels, "target")

    if target_profile_error
        if target_profile_model == :dwarf_Zhao
            params = product(
                [target_profile_ρ0, target_profile_ρ0_u, target_profile_ρ0_d],
                [target_profile_rs, target_profile_rs_u, target_profile_rs_d],
                [target_profile_α, target_profile_α_u, target_profile_α_d],
                [target_profile_β, target_profile_β_u, target_profile_β_d],
                [target_profile_γ, target_profile_γ_u, target_profile_γ_d],
            )
            models = [Zhao(p...) for p in params]
            ρ_halo_target_u = similar(ρ_halo_target)
            ρ_halo_target_d = similar(ρ_halo_target)
            for i in eachindex(r_target)
                rhos = ustrip.(u"Msun/kpc^3", GalacticDynamics.density.(models, r_target[i] * u"kpc"))
                u, d = extrema(rhos)
                ρ_halo_target_u[i] = u
                ρ_halo_target_d[i] = d
            end
        elseif target_profile_model == :dwarf_gNFW
            params = product(
                [target_profile_β, target_profile_β_u, target_profile_β_d],
                [target_profile_ρ0, target_profile_ρ0_u, target_profile_ρ0_d],
                [target_profile_rs, target_profile_rs_u, target_profile_rs_d],
            )
            models = [gNFW(p...) for p in params]
            ρ_halo_target_u = similar(ρ_halo_target)
            ρ_halo_target_d = similar(ρ_halo_target)
            for i in eachindex(r_target)
                rhos = ustrip.(u"Msun/kpc^3", GalacticDynamics.density.(models, r_target[i] * u"kpc"))
                u, d = extrema(rhos)
                ρ_halo_target_u[i] = u
                ρ_halo_target_d[i] = d
            end
        end
        lt_u = Makie.lines!(AxisDensityProfile, r_target, ρ_halo_target_u, color = (:black, 0.3))
        lt_d = Makie.lines!(AxisDensityProfile, r_target, ρ_halo_target_d, color = (:black, 0.3))
    end

    Legend(fig[2,3], plots, labels;
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
        margin = (10, 10, 10, 10),
    )
    
    return (
        profile_r_mean = profile_r_mean,
        profile_ρ_mean = profile_ρ_mean,
        r_target = r_target,
        ρ_halo_target = ρ_halo_target,
    )
end

"""
$(TYPEDSIGNATURES)

Update progress bar with Unicode plots if enabled.
"""
function update_unicode_progress!(progress, i, t, unicode_plot, distributed_memory, rho, rho_max_id,
    Realtime, StepsBetweenSnapshots, r_target, ρ_halo_target, _profile_r_mean, _profile_ρ_mean,
    best_fit_t, best_fit_error, current_fit_error, best_fit_beta_star_error, best_fit_beta_star, current_beta_star,
    unicode_heatmap_width, Xmax, uT, uL, Nx, Δ)
    
    if unicode_plot
        if distributed_memory
            slice_rho = dropdims(collect(rho[:, :, rho_max_id[3]]), dims=3)
        else
            slice_rho = rho[:, :, rho_max_id[3]]
        end

        unicode_p1 = UnicodePlots.heatmap(log10.(slice_rho);
            xoffset = -Xmax*uL,
            yoffset = -Ymax*uL,
            xfact = 2*Xmax*uL/Nx,
            yfact = 2*Xmax*uL/Nx,
            height = unicode_heatmap_width,
            width = unicode_heatmap_width,
        )

        if Realtime && iszero(mod(i, StepsBetweenSnapshots))
            unicode_p2 = UnicodePlots.lineplot(r_target, ρ_halo_target;
                xscale=:log10,
                yscale=:log10,
                xlim = (0.5*Δ[1]*uL, Xmax*uL),
                ylim = (1e3, 1e10),
                color = :blue,
                xlabel = "log10(r [kpc])",
                ylabel = "log10(ρ [Msun/kpc³])",
            )
            UnicodePlots.lineplot!(unicode_p2, _profile_r_mean, _profile_ρ_mean;
                color = :red,
            )
        else
            unicode_p2 = nothing
        end

        next!(progress; showvalues = [
            ("iter", i),
            ("t [Gyr]", t[i] * uT),
            ("density", unicode_p1),
            ("profile", unicode_p2),
            ("best_fit_t", best_fit_t),
            ("best_fit_error", best_fit_error),
            ("current_fit_error", current_fit_error),
            ("best_fit_beta_star_error", best_fit_beta_star_error),
            ("best_fit_beta_star", best_fit_beta_star),
            ("current_beta_star", current_beta_star),
        ])
    else
        next!(progress; showvalues = [
            ("iter", i),
            ("t [Gyr]", t[i] * uT),
            ("best_fit_t", best_fit_t),
            ("best_fit_error", best_fit_error),
            ("current_fit_error", current_fit_error),
            ("best_fit_beta_star_error", best_fit_beta_star_error),
            ("best_fit_beta_star", best_fit_beta_star),
            ("current_beta_star", current_beta_star),
        ])
    end
end

"""
$(TYPEDSIGNATURES)

Plot MOND acceleration comparison.
This function encapsulates the `plotMOND` inner function from SPE3D_MOND.
"""
function plotMOND(ax_all, ay_all, az_all, ax_b, ay_b, az_b, a0, r, length_astro, acc_astro, minR, maxR, outputdir, title, suffix, section;
    filename = title)
    
    rMOND = r[:, :, div(end,2)][:] * ustrip(length_astro)
    a_b = sqrt.(ax_b[:, :, div(end,2)].^2 .+ ay_b[:, :, div(end,2)].^2 .+ az_b[:, :, div(end,2)].^2)
    a_all = sqrt.(ax_all[:, :, div(end,2)].^2 .+ ay_all[:, :, div(end,2)].^2 .+ az_all[:, :, div(end,2)].^2)

    a_mond = a_b ./ (1 .- exp.(-sqrt.(a_b./a0)))

    figMOND = Figure(size = (1600, 900));
    ax = Axis(figMOND[1,1];
        xlabel = "r [kpc]",
        ylabel = "acc [m/s²]",
    )

    r_mean, a_b_mean, r_std, a_b_std = distribution(rMOND[:], a_b[:]; section)
    _, a_all_mean, _, a_all_std = distribution(rMOND[:], a_all[:];    section)
    _, a_mond_mean, _, a_mond_std = distribution(rMOND[:], a_mond[:]; section)

    # relative error of accelerations in region [minR, maxR]
    indexMinR = findfirst(x->x>ustrip(u"kpc", minR), r_mean)
    indexMaxR = findfirst(x->x>ustrip(u"kpc", maxR), r_mean)
    a_all_measurement = measurement.(a_all_mean[indexMinR:indexMaxR], a_all_std[indexMinR:indexMaxR])
    a_mond_measurement = measurement.(a_mond_mean[indexMinR:indexMaxR], a_mond_std[indexMinR:indexMaxR])
    MOND_errorrel = mean((a_all_measurement .- a_mond_measurement) ./ a_mond_measurement)
    
    uAcc = ustrip(acc_astro)

    f1 = Makie.scatter!(ax, rMOND[:], a_b[:] * uAcc, markersize = 2, color = :red)
    f1a = Makie.lines!(ax, r_mean, a_b_mean * uAcc, color = :red)
    f1b = Makie.band!(ax, r_mean, (a_b_mean - a_b_std) * uAcc, (a_b_mean + a_b_std) * uAcc, color = (:red, 0.2))

    f2 = Makie.scatter!(ax, rMOND[:], a_all[:] * uAcc, markersize = 2, color = :black)
    f2a = Makie.lines!(ax, r_mean, a_all_mean * uAcc, color = :black)
    f2b = Makie.band!(ax, r_mean, (a_all_mean - a_all_std) * uAcc, (a_all_mean + a_all_std) * uAcc, color = (:black, 0.2))

    f3 = Makie.scatter!(ax, rMOND[:], a_mond[:] * uAcc, markersize = 2, color = :green)
    f3a = Makie.lines!(ax, r_mean, a_mond_mean * uAcc, color = :green)
    f3b = Makie.band!(ax, r_mean, (a_mond_mean - a_mond_std) * uAcc, (a_mond_mean + a_mond_std) * uAcc, color = (:green, 0.2))

    f4 = Makie.hlines!(ax, [a0 * uAcc], color = :blue)
    
    Legend(figMOND[1,2], [[f1, f1a, f1b], [f2, f2a, f2b], [f3, f3a, f3b], f4], ["baryon", "all", "mond", "a0"])
    Makie.xlims!(ax, 0, ustrip(u"kpc", 1.5 * maxR))
    Makie.ylims!(ax, 1.0e-11, 1.0e-9)
    Makie.save(joinpath(outputdir, filename * "- RAR, $(suffix).png"), figMOND)
    return figMOND, MOND_errorrel
end

# Export functions
export setup_visualization, setup_virial_visualization
export update_virial_terms!, update_virial_visualization!, setup_density_profile_visualization!
export update_unicode_progress!, plotMOND
