# SPE visualization module

"""
$(TYPEDSIGNATURES)

Setup initial visualization figure with multiple panels.
"""
function setup_visualization(Xmax, t, Δ, rho, rho_max_id, total_halo_mass, radii,
    uT, uL, title, suffix, distributed_memory;
    size = (2400, 1400)
)
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
    ArrayR = Observable([std(collect(r[:,:,div(end,2)]), aweights(collect(rho[:,:,div(end,2)]))) * uL]) #TODO this is strongly affected by rho
    ArrayR1 = Observable([radii[1] * uL])
    ArrayR2 = Observable([radii[2] * uL])
    ArrayR3 = Observable([radii[3] * uL])
    ArrayR4 = Observable([radii[4] * uL])
    ArrayR5 = Observable([radii[5] * uL])
    ArrayR6 = Observable([radii[6] * uL])
    ArrayR7 = Observable([radii[7] * uL])
    ArrayR8 = Observable([radii[8] * uL])
    ArrayR9 = Observable([radii[9] * uL])

    ColorRange = Observable((0, maximum(abs.(ψ).^2)/10))
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
        fig = fig,
        ArrayT = ArrayT,
        ArrayT_Snap = ArrayT_Snap,
        AxisR = AxisR,
        AxisVirial = AxisVirial,
        AxisDensityProfile = AxisDensityProfile,
        SliceXY = SliceXY,
        SliceYZ = SliceYZ,
        SliceXZ = SliceXZ,
        ArrayTotalMass = ArrayTotalMass,
        ArrayR = ArrayR,
        ArrayR1 = ArrayR1,
        ArrayR2 = ArrayR2,
        ArrayR3 = ArrayR3,
        ArrayR4 = ArrayR4,
        ArrayR5 = ArrayR5,
        ArrayR6 = ArrayR6,
        ArrayR7 = ArrayR7,
        ArrayR8 = ArrayR8,
        ArrayR9 = ArrayR9,
        ColorRange = ColorRange,
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

# Export functions
export setup_visualization, setup_virial_visualization, update_visualization!
export update_virial_terms!, update_virial_visualization!, setup_density_profile_visualization!
