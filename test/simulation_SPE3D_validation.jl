# ============================================================================
# 3D Schrödinger-Poisson Equation Validation Tests
# Using WaveDM.jl's SPE3D_waveDM function
# Based on PyUltraLight (Edwards et al. 2018)
# ============================================================================
#
# Usage: include("test/simulation_SPE3D_validation.jl")
#
# This will automatically run all convergence & validation tests and generate
# publication-quality figures saved to output/SPE3D/figures/
#
# ============================================================================

using WaveDM
using LinearAlgebra
using FFTW
using ProgressMeter
using Statistics
using Printf
using CairoMakie
using FileIO

const OUTPUT_DIR = "./output/SPE3D"
const FIGURES_DIR = joinpath(OUTPUT_DIR, "figures")

mkpath.([OUTPUT_DIR, FIGURES_DIR])

# ============================================================================
# Analytic Solutions
# ============================================================================

function analytic_soliton_3d(x, y, z, t; A=sqrt(2.0), v=[0.0, 0.0, 0.0], x0=[0.0, 0.0, 0.0])
    xi_x = x .- x0[1] .- v[1] .* t
    xi_y = y .- x0[2] .- v[2] .* t
    xi_z = z .- x0[3] .- v[3] .* t
    r = sqrt.(xi_x.^2 .+ xi_y.^2 .+ xi_z.^2)
    B = A / 2
    omega = A^2 / 2 - sum(v.^2) / 2
    phase = v[1] .* xi_x .+ v[2] .* xi_y .+ v[3] .* xi_z .- omega .* t
    return A .* sech.(B .* r) .* exp.(im .* phase)
end

function analytic_soliton_collision_3d(x, y, z, t; 
    A1=sqrt(2.0), v1=[0.3, 0.0, 0.0], A2=sqrt(2.0), v2=[-0.3, 0.0, 0.0], 
    x1=[15.0, 0.0, 0.0], x2=[-15.0, 0.0, 0.0])
    return analytic_soliton_3d(x, y, z, t; A=A1, v=v1, x0=x1) .+ 
           analytic_soliton_3d(x, y, z, t; A=A2, v=v2, x0=x2)
end

function analytic_soliton_binary_3d(x, y, z, t; 
    A=sqrt(2.0), separation=20.0, orbital_velocity=0.2)
    omega_orb = orbital_velocity / (separation / 2)
    x1 = [separation/2 * cos(omega_orb * t), separation/2 * sin(omega_orb * t), 0.0]
    x2 = [-separation/2 * cos(omega_orb * t), -separation/2 * sin(omega_orb * t), 0.0]
    v1 = [-orbital_velocity * sin(omega_orb * t), orbital_velocity * cos(omega_orb * t), 0.0]
    v2 = [orbital_velocity * sin(omega_orb * t), -orbital_velocity * cos(omega_orb * t), 0.0]
    return analytic_soliton_3d(x, y, z, t; A=A, v=v1, x0=x1) .+ 
           analytic_soliton_3d(x, y, z, t; A=A, v=v2, x0=x2)
end

# ============================================================================
# Error Metrics
# ============================================================================

compute_error_l2(psi_num, psi_ana) = sqrt(mean(abs.(psi_num .- psi_ana).^2))
compute_error_linf(psi_num, psi_ana) = maximum(abs.(psi_num .- psi_ana))
compute_relative_error_l2(psi_num, psi_ana) =
    sqrt(mean(abs.(psi_num .- psi_ana).^2)) / sqrt(mean(abs.(psi_ana).^2))

# ============================================================================
# Energy Computation (simplified for validation - no internal functions)
# ============================================================================

function compute_total_mass(ψ, Δ)
    return sum(abs2.(ψ)) * prod(Δ)
end

function compute_peak_density(ψ)
    return maximum(abs2.(ψ))
end

function find_soliton_center(ψ, x, y, z)
    rho = abs2.(ψ)
    idx = argmax(rho)
    ix, iy, iz = Tuple(idx)
    cx = x[ix]
    cy = y[iy]
    cz = z[iz]
    return cx, cy, cz
end

function compute_radial_profile(ψ, x, y, z)
    Nx = length(x)
    iz_center = div(Nx, 2)
    iy_center = div(Nx, 2)
    rho_slice = abs2.(ψ[:, iy_center, :])
    r = sqrt.(x.^2 .+ z.^2')
    return r, rho_slice
end

# ============================================================================
# Test Functions
# ============================================================================

function run_single_soliton_test(; Xmax=30.0, Tmax=2.0, Nx=64, Nt=64, 
    A=sqrt(2.0), v=[0.1, 0.0, 0.0], x0=[0.0, 0.0, 0.0], verbose=true)
    
    verbose && @info "Running single soliton test: Nx=$Nx, Nt=$Nt"
    
    IC = (x, y, z) -> analytic_soliton_3d(x, y, z, 0.0; A=A, v=v, x0=x0)
    
    ψ_final, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2 = SPE3D_waveDM(;
        Xmax=Xmax, Ymax=Xmax, Zmax=Xmax, Tmax=Tmax,
        Nx=Nx, Ny=Nx, Nz=Nx, Nt=Nt,
        IC=IC,
        V=(x, y, z, ψ) -> 0.0,
        baryon_mode=:ignored,
        absorb_coeff=0.0,
        outputdir=joinpath(OUTPUT_DIR, "single_soliton"),
        title="Single Soliton",
        gpu=false,
        Realtime=false,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
    )

    x = collect(LinRange(-Xmax, Xmax, Nx))
    y = collect(LinRange(-Xmax, Xmax, Nx))
    z = collect(LinRange(-Xmax, Xmax, Nx))
    xxx = [x[i] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    yyy = [y[j] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    zzz = [z[k] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    
    psi_ana = analytic_soliton_3d(xxx, yyy, zzz, Tmax; A=A, v=v, x0=x0)
    
    return Dict(
        :ψ => ψ_final,
        :x => x, :y => y, :z => z,
        :Nx => Nx, :Nt => Nt, :Xmax => Xmax, :Tmax => Tmax,
        :A => A, :v => v, :x0 => x0,
        :psi_ana => psi_ana,
        :error_l2 => compute_error_l2(ψ_final, psi_ana),
        :error_linf => compute_error_linf(ψ_final, psi_ana),
        :rel_error_l2 => compute_relative_error_l2(ψ_final, psi_ana)
    )
end

function run_soliton_collision_test(; Xmax=60.0, Tmax=10.0, Nx=64, Nt=64, verbose=true)
    
    verbose && @info "Running soliton collision test: Nx=$Nx, Nt=$Nt"
    
    IC = (x, y, z) -> analytic_soliton_collision_3d(x, y, z, 0.0)
    
    ψ_final, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2 = SPE3D_waveDM(;
        Xmax=Xmax, Ymax=Xmax, Zmax=Xmax, Tmax=Tmax,
        Nx=Nx, Ny=Nx, Nz=Nx, Nt=Nt,
        IC=IC,
        V=(x, y, z, ψ) -> 0.0,
        baryon_mode=:ignored,
        absorb_coeff=0.0,
        outputdir=joinpath(OUTPUT_DIR, "soliton_collision"),
        title="Soliton Collision",
        gpu=false,
        Realtime=false,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
    )

    return Dict(
        :ψ => ψ_final,
        :x => collect(LinRange(-Xmax, Xmax, Nx)),
        :y => collect(LinRange(-Xmax, Xmax, Nx)),
        :z => collect(LinRange(-Xmax, Xmax, Nx)),
        :Nx => Nx, :Nt => Nt, :Xmax => Xmax, :Tmax => Tmax
    )
end

function analytic_soliton_tde_3d(x, y, z, t; A=sqrt(2.0), x0=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0])
    return analytic_soliton_3d(x, y, z, t; A=A, v=v, x0=x0)
end

function run_tde_test(; Xmax=50.0, Tmax=15.0, Nx=64, Nt=64,
    A=sqrt(2.0), soliton_x0=[0.0, 0.0, 0.0], soliton_v=[0.0, 0.0, 0.0],
    perturber_mass=0.5, perturber_x0=[25.0, 0.0, 0.0], perturber_v=[-0.2, 0.0, 0.0],
    verbose=true)

    verbose && @info "Running TDE test: M_perturber=$perturber_mass"

    function IC_with_perturber(x, y, z)
        soliton = analytic_soliton_3d(x, y, z, 0.0; A=A, v=soliton_v, x0=soliton_x0)
        r_perturber = sqrt.((x .- perturber_x0[1]).^2 + (y .- perturber_x0[2]).^2 + (z .- perturber_x0[3]).^2)
        gaussian_perturber = perturber_mass * exp.(-r_perturber.^2 / (2 * 3.0^2))
        return soliton + gaussian_perturber
    end

    function V_with_gravity(x, y, z, ψ)
        r_perturber = sqrt.((x .- perturber_x0[1]).^2 + (y .- perturber_x0[2]).^2 + (z .- perturber_x0[3]).^2)
        G = 1.0
        M = perturber_mass
        softening = 3.0
        potential = -G * M ./ sqrt.(r_perturber.^2 .+ softening^2)
        return potential
    end

    ψ_final, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2 = SPE3D_waveDM(;
        Xmax=Xmax, Ymax=Xmax, Zmax=Xmax, Tmax=Tmax,
        Nx=Nx, Ny=Nx, Nz=Nx, Nt=Nt,
        IC=IC_with_perturber,
        V=V_with_gravity,
        baryon_mode=:ignored,
        absorb_coeff=0.0,
        outputdir=joinpath(OUTPUT_DIR, "tde"),
        title="TDE",
        gpu=false,
        Realtime=false,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
    )

    x = collect(LinRange(-Xmax, Xmax, Nx))
    y = collect(LinRange(-Xmax, Xmax, Nx))
    z = collect(LinRange(-Xmax, Xmax, Nx))
    xxx = [x[i] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    yyy = [y[j] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    zzz = [z[k] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    Δ = ntuple(_ -> 2*Xmax/Nx, 3)

    ψ_initial = IC_with_perturber(xxx, yyy, zzz)
    mass_initial = compute_total_mass(ψ_initial, Δ)
    mass_final = compute_total_mass(ψ_final, Δ)

    return Dict(
        :ψ => ψ_final,
        :x => x, :y => y, :z => z,
        :Nx => Nx, :Nt => Nt, :Xmax => Xmax, :Tmax => Tmax,
        :mass_ratio => mass_final / mass_initial,
        :mass_initial => mass_initial,
        :mass_final => mass_final,
        :perturber_mass => perturber_mass
    )
end

function run_soliton_binary_test(; Xmax=50.0, Tmax=20.0, Nx=64, Nt=64, verbose=true)
    
    verbose && @info "Running soliton binary test: Nx=$Nx, Nt=$Nt"
    
    IC = (x, y, z) -> analytic_soliton_binary_3d(x, y, z, 0.0)
    
    ψ_final, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2 = SPE3D_waveDM(;
        Xmax=Xmax, Ymax=Xmax, Zmax=Xmax, Tmax=Tmax,
        Nx=Nx, Ny=Nx, Nz=Nx, Nt=Nt,
        IC=IC,
        V=(x, y, z, ψ) -> 0.0,
        baryon_mode=:ignored,
        absorb_coeff=0.0,
        outputdir=joinpath(OUTPUT_DIR, "soliton_binary"),
        title="Soliton Binary",
        gpu=false,
        Realtime=false,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
    )

    return Dict(
        :ψ => ψ_final,
        :x => collect(LinRange(-Xmax, Xmax, Nx)),
        :y => collect(LinRange(-Xmax, Xmax, Nx)),
        :z => collect(LinRange(-Xmax, Xmax, Nx)),
        :Nx => Nx, :Nt => Nt, :Xmax => Xmax, :Tmax => Tmax
    )
end

function run_tde_convergence_test(; Xmax=50.0, Tmax=15.0, Nx_list=[48, 64, 96, 128], Nt_list=nothing, verbose=true)
    if isnothing(Nt_list)
        Nt_list = Nx_list
    end

    verbose && @info "Running TDE convergence test"

    profiles = []
    mass_ratios = Float64[]

    for (i, Nx) in enumerate(Nx_list)
        verbose && @info "  Nx = $Nx"
        Nt = Nt_list[min(i, length(Nt_list))]
        result = run_tde_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        push!(mass_ratios, result[:mass_ratio])

        x = result[:x]
        N = length(x)
        iz_center = div(N, 2)
        iy_center = div(N, 2)
        rho_line = abs2.(result[:ψ])[iz_center, iy_center, :]
        push!(profiles, (x=result[:x], rho=rho_line, label="Nx=$Nx"))
    end

    return Dict(
        :Nx_list => Nx_list,
        :profiles => profiles,
        :mass_ratios => mass_ratios
    )
end

function run_collision_convergence_test(; Xmax=60.0, Tmax=10.0, Nx_list=[48, 64, 96, 128], Nt_list=nothing, verbose=true)
    if isnothing(Nt_list)
        Nt_list = Nx_list
    end

    verbose && @info "Running collision convergence test"

    profiles = []

    for (i, Nx) in enumerate(Nx_list)
        verbose && @info "  Nx = $Nx"
        Nt = Nt_list[min(i, length(Nt_list))]
        result = run_soliton_collision_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        x = result[:x]
        N = length(x)
        iz_center = div(N, 2)
        rho_line = abs2.(result[:ψ])[:, :, iz_center]
        push!(profiles, (x=result[:x], y=result[:y], rho=rho_line, label="Nx=$Nx"))
    end

    return Dict(
        :Nx_list => Nx_list,
        :profiles => profiles
    )
end

function run_binary_convergence_test(; Xmax=50.0, Tmax=20.0, Nx_list=[48, 64, 96, 128], Nt_list=nothing, verbose=true)
    if isnothing(Nt_list)
        Nt_list = Nx_list
    end

    verbose && @info "Running binary convergence test"

    profiles = []

    for (i, Nx) in enumerate(Nx_list)
        verbose && @info "  Nx = $Nx"
        Nt = Nt_list[min(i, length(Nt_list))]
        result = run_soliton_binary_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        x = result[:x]
        N = length(x)
        iz_center = div(N, 2)
        rho_line = abs2.(result[:ψ])[:, :, iz_center]
        push!(profiles, (x=result[:x], y=result[:y], rho=rho_line, label="Nx=$Nx"))
    end

    return Dict(
        :Nx_list => Nx_list,
        :profiles => profiles
    )
end

function run_spatial_convergence_test(; Xmax=30.0, Tmax=2.0, Nx_list=[32, 64, 128, 192], Nt=64, verbose=true)
    
    verbose && @info "Running spatial convergence test"
    
    errors_l2 = Float64[]
    errors_linf = Float64[]
    rel_errors_l2 = Float64[]
    profiles = []

    for Nx in Nx_list
        verbose && @info "  Nx = $Nx"
        result = run_single_soliton_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)
        push!(errors_l2, result[:error_l2])
        push!(errors_linf, result[:error_linf])
        push!(rel_errors_l2, result[:rel_error_l2])
        push!(profiles, (x=result[:x], rho=abs2.(result[:ψ]), label="Nx=$Nx"))
    end
    
    dx_list = 2 * Xmax ./ Nx_list
    orders_l2 = Float64[]
    orders_linf = Float64[]
    
    for i in 2:length(Nx_list)
        if errors_l2[i-1] > 0
            push!(orders_l2, log(errors_l2[i]/errors_l2[i-1]) / log(dx_list[i]/dx_list[i-1]))
        else
            push!(orders_l2, NaN)
        end
        if errors_linf[i-1] > 0
            push!(orders_linf, log(errors_linf[i]/errors_linf[i-1]) / log(dx_list[i]/dx_list[i-1]))
        else
            push!(orders_linf, NaN)
        end
    end
    
    return Dict(
        :Nx_list => Nx_list,
        :dx_list => dx_list,
        :errors_l2 => errors_l2,
        :errors_linf => errors_linf,
        :rel_errors_l2 => rel_errors_l2,
        :orders_l2 => orders_l2,
        :orders_linf => orders_linf,
        :profiles => profiles
    )
end

function run_temporal_convergence_test(; Xmax=30.0, Tmax=2.0, Nx=64, Nt_list=[32, 64, 128, 192], verbose=true)

    verbose && @info "Running temporal convergence test"

    errors_l2 = Float64[]
    errors_linf = Float64[]
    rel_errors_l2 = Float64[]
    profiles = []

    for Nt in Nt_list
        verbose && @info "  Nt = $Nt"
        result = run_single_soliton_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)
        push!(errors_l2, result[:error_l2])
        push!(errors_linf, result[:error_linf])
        push!(rel_errors_l2, result[:rel_error_l2])
        push!(profiles, (x=result[:x], rho=abs2.(result[:ψ]), label="Nt=$Nt"))
    end
    
    dt_list = Tmax ./ Nt_list
    orders_l2 = Float64[]
    orders_linf = Float64[]
    
    for i in 2:length(Nt_list)
        if errors_l2[i-1] > 0
            push!(orders_l2, log(errors_l2[i]/errors_l2[i-1]) / log(dt_list[i]/dt_list[i-1]))
        else
            push!(orders_l2, NaN)
        end
        if errors_linf[i-1] > 0
            push!(orders_linf, log(errors_linf[i]/errors_linf[i-1]) / log(dt_list[i]/dt_list[i-1]))
        else
            push!(orders_linf, NaN)
        end
    end
    
    return Dict(
        :Nt_list => Nt_list,
        :dt_list => dt_list,
        :errors_l2 => errors_l2,
        :errors_linf => errors_linf,
        :rel_errors_l2 => rel_errors_l2,
        :orders_l2 => orders_l2,
        :orders_linf => orders_linf,
        :profiles => profiles
    )
end

# ============================================================================
# Publication-Quality Plotting Functions
# ============================================================================

function setup_font!(ax; xlabel="", ylabel="", title="", titlesize=14, labelsize=12, ticksize=10)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.title = title
    ax.titlesize = titlesize
    ax.labelsize = labelsize
    ax.ticksize = ticksize
    ax.xtickwidth = 1.5
    ax.ytickwidth = 1.5
    ax.spineswidth = 1.5
end

function plot_spatial_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_spatial_convergence.png"))

    Nx_list = result[:Nx_list]
    dx_list = result[:dx_list]
    errors_l2 = result[:errors_l2]
    errors_linf = result[:errors_linf]
    orders_l2 = result[:orders_l2]
    profiles = result[:profiles]

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: Density profiles overlay
    ax1 = Axis(fig[1, 1],
        title="Density Profile Convergence",
        xlabel="r",
        ylabel="|ψ|²",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red, :purple]
    alphas = [0.7, 0.8, 0.9, 1.0, 1.0]
    for (i, p) in enumerate(profiles)
        Nx = length(p.x)
        iz_center = div(Nx, 2)
        iy_center = div(Nx, 2)
        rho_line = p.rho[iz_center, iy_center, :]
        r_line = abs.(p.x)
        lines!(ax1, r_line, rho_line, color=colors[i], linewidth=2, alpha=alphas[i], label=p.label)
    end
    axislegend(ax1, position=:rt, framevisible=true)
    xlims!(ax1, 0, 15)

    # Right panel: Error convergence
    ax2 = Axis(fig[1, 2],
        title="Spatial Convergence Order",
        xlabel="Grid spacing h",
        ylabel="L2 Error",
        xscale=log10,
        yscale=log10,
        xgridvisible=true,
        ygridvisible=true)

    scatter!(ax2, dx_list, errors_l2,
        color=:blue, marker=:circle, markersize=10,
        label="Numerical error")

    # Reference lines
    h_ref = dx_list[1]
    order2_line = errors_l2[1] * (dx_list ./ h_ref).^2
    order4_line = errors_l2[1] * (dx_list ./ h_ref).^4

    lines!(ax2, dx_list, order2_line,
        color=:gray, linestyle=:dash, linewidth=1.5, label="O(h²)")
    lines!(ax2, dx_list, order4_line,
        color=:lightgray, linestyle=:dot, linewidth=1.5, label="O(h⁴)")

    axislegend(ax2, position=:lb, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    println("\n" * "="^60)
    println("Spatial Convergence Results:")
    println("="^60)
    println("Nx\t\th\t\tL2 Error\t\tL∞ Error\t\tOrder L2")
    for i in eachindex(Nx_list)
        order_str = i > 1 ? @sprintf("%.2f", orders_l2[i-1]) : "-"
        println("$(Nx_list[i])\t\t$(@sprintf("%.4e", dx_list[i]))\t\t$(@sprintf("%.4e", errors_l2[i]))\t\t$(@sprintf("%.4e", errors_linf[i]))\t\t$order_str")
    end
    println("="^60)

    return fig
end

function plot_temporal_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_temporal_convergence.png"))

    Nt_list = result[:Nt_list]
    dt_list = result[:dt_list]
    errors_l2 = result[:errors_l2]
    errors_linf = result[:errors_linf]
    orders_l2 = result[:orders_l2]
    profiles = result[:profiles]

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: Density profiles overlay
    ax1 = Axis(fig[1, 1],
        title="Temporal Convergence - Density Profiles",
        xlabel="r",
        ylabel="|ψ|²",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red, :purple]
    alphas = [0.7, 0.8, 0.9, 1.0, 1.0]
    for (i, p) in enumerate(profiles)
        Nx = length(p.x)
        iz_center = div(Nx, 2)
        iy_center = div(Nx, 2)
        rho_line = p.rho[iz_center, iy_center, :]
        r_line = abs.(p.x)
        lines!(ax1, r_line, rho_line, color=colors[i], linewidth=2, alpha=alphas[i], label=p.label)
    end
    axislegend(ax1, position=:rt, framevisible=true)
    xlims!(ax1, 0, 15)

    # Right panel: Error convergence
    ax2 = Axis(fig[1, 2],
        title="Temporal Convergence Order",
        xlabel="Time step Δt",
        ylabel="L2 Error",
        xscale=log10,
        yscale=log10,
        xgridvisible=true,
        ygridvisible=true)

    scatter!(ax2, dt_list, errors_l2,
        color=:blue, marker=:circle, markersize=10,
        label="Numerical error")

    # Reference lines
    dt_ref = dt_list[1]
    order2_line = errors_l2[1] * (dt_list ./ dt_ref).^2

    lines!(ax2, dt_list, order2_line,
        color=:gray, linestyle=:dash, linewidth=1.5, label="O(Δt²)")

    axislegend(ax2, position=:lb, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    println("\n" * "="^60)
    println("Temporal Convergence Results:")
    println("="^60)
    println("Nt\t\tΔt\t\tL2 Error\t\tL∞ Error\t\tOrder L2")
    for i in eachindex(Nt_list)
        order_str = i > 1 ? @sprintf("%.2f", orders_l2[i-1]) : "-"
        println("$(Nt_list[i])\t\t$(@sprintf("%.4e", dt_list[i]))\t\t$(@sprintf("%.4e", errors_l2[i]))\t\t$(@sprintf("%.4e", errors_linf[i]))\t\t$order_str")
    end
    println("="^60)

    return fig
end

function plot_tde_result(result; save_path=joinpath(FIGURES_DIR, "fig_tde.png"))

    ψ = result[:ψ]
    x = result[:x]
    y = result[:y]
    z = result[:z]
    Nx = result[:Nx]
    rho = abs2.(ψ)
    rho_max = maximum(rho)
    iz_center = div(Nx, 2)

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: XY slice
    ax1 = Axis(fig[1, 1],
        title="TDE - XY Slice (z=0)",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm1 = heatmap!(ax1, x, y, rho[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 2], hm1)

    # Right panel: XZ slice
    ax2 = Axis(fig[1, 3],
        title="TDE - XZ Slice (y=0)",
        xlabel="x",
        ylabel="z",
        xgridvisible=true,
        ygridvisible=true)
    hm2 = heatmap!(ax2, x, z, rho[:, iz_center, :],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 4], hm2)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    println("\n" * "="^60)
    println("TDE Test Results:")
    println("="^60)
    println("Initial Mass: $(@sprintf("%.6f", result[:mass_initial]))")
    println("Final Mass:   $(@sprintf("%.6f", result[:mass_final]))")
    println("Mass Ratio:   $(@sprintf("%.2f", 100*result[:mass_ratio]))%")
    println("Perturber Mass: $(result[:perturber_mass])")
    println("="^60)

    return fig
end

function plot_energy_conservation(result; save_path=joinpath(FIGURES_DIR, "fig_energy_conservation.png"))

    ψ = result[:ψ]
    Nx = result[:Nx]
    Xmax = result[:Xmax]

    Δ = 2 * Xmax / Nx
    mass = compute_total_mass(ψ, ntuple(_ -> Δ, 3))
    peak_density = compute_peak_density(ψ)
    cx, cy, cz = find_soliton_center(ψ, result[:x], result[:y], result[:z])

    println("\n" * "="^60)
    println("Soliton Propagation Metrics:")
    println("="^60)
    println("Total Mass: $(@sprintf("%.6f", mass))")
    println("Peak Density: $(@sprintf("%.6f", peak_density))")
    println("Soliton Center: ($(@sprintf("%.3f", cx)), $(@sprintf("%.3f", cy)), $(@sprintf("%.3f", cz)))")
    println("Expected Center: (0.0, 0.0, 0.0)")
    println("Center Error: $(@sprintf("%.6f", sqrt(cx^2 + cy^2 + cz^2)))")
    println("="^60)

    fig = Figure(size=(600, 400), fontsize=12)
    ax = Axis(fig[1, 1], title="Soliton Propagation Metrics")

    # Create a simple table-like display
    metrics = [
        "Total Mass" => @sprintf("%.4f", mass),
        "Peak Density" => @sprintf("%.4f", peak_density),
        "Center X" => @sprintf("%.3f", cx),
        "Center Y" => @sprintf("%.3f", cy),
        "Center Z" => @sprintf("%.3f", cz),
        "Center Error" => @sprintf("%.4f", sqrt(cx^2+cy^2+cz^2))
    ]

    y_pos = 0.9
    for (label, value) in metrics
        text!(ax, Point2(0.1, y_pos), text="$label:", fontsize=12, align=(:left, :top))
        text!(ax, Point2(0.6, y_pos), text=value, fontsize=12, align=(:left, :top))
        y_pos -= 0.12
    end

    hidedecorations!(ax)
    hidespines!(ax)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_soliton_density(result; save_path=joinpath(FIGURES_DIR, "fig_single_soliton_density.png"))

    ψ = result[:ψ]
    x = result[:x]
    y = result[:y]
    z = result[:z]
    Nx = result[:Nx]
    rho = abs2.(ψ)
    rho_max = maximum(rho)
    iz_center = div(Nx, 2)

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: XY slice
    ax1 = Axis(fig[1, 1],
        title="Soliton Density - XY Slice",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm1 = heatmap!(ax1, x, y, rho[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 2], hm1)

    # Right panel: Profile comparison
    ax2 = Axis(fig[1, 3],
        title="Radial Profile Comparison",
        xlabel="r",
        ylabel="|ψ|²",
        xgridvisible=true,
        ygridvisible=true)

    x_line = x
    rho_num_line = rho[iz_center, iz_center, :]
    r_line = abs.(x_line)
    rho_ana_line = abs2.(analytic_soliton_3d(x_line, 0, 0, 0; A=result[:A], v=result[:v], x0=result[:x0]))

    lines!(ax2, r_line, rho_num_line, color=:blue, linewidth=2.5, label="Numerical")
    lines!(ax2, r_line, rho_ana_line, color=:red, linewidth=2, linestyle=:dash, label="Analytic")
    axislegend(ax2, position=:rt, framevisible=true)
    xlims!(ax2, 0, 15)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    println("\n" * "="^60)
    println("Single Soliton Test Results:")
    println("="^60)
    println("Nx = $(result[:Nx]), Nt = $(result[:Nt])")
    println("L2 Error:        $(@sprintf("%.6e", result[:error_l2]))")
    println("L∞ Error:        $(@sprintf("%.6e", result[:error_linf]))")
    println("Relative L2 Error: $(@sprintf("%.6e", result[:rel_error_l2]))")
    println("="^60)

    return fig
end

function plot_soliton_collision(result; save_path=joinpath(FIGURES_DIR, "fig_soliton_collision.png"))

    ψ = result[:ψ]
    x = result[:x]
    y = result[:y]
    z = result[:z]
    Nx = result[:Nx]
    Tmax = result[:Tmax]
    rho_final = abs2.(ψ)
    rho_max = maximum(rho_final)
    iz_center = div(Nx, 2)

    xxx = [x[i] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    yyy = [y[j] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    zzz = [z[k] for i in 1:Nx, j in 1:Nx, k in 1:Nx]

    rho_t0 = abs2.(analytic_soliton_collision_3d(xxx, yyy, zzz, 0.0))
    rho_tmid = abs2.(analytic_soliton_collision_3d(xxx, yyy, zzz, Tmax/2))

    fig = Figure(size=(1200, 400), fontsize=11)

    # Initial state
    ax1 = Axis(fig[1, 1],
        title="Initial State (t=0)",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm1 = heatmap!(ax1, x, y, rho_t0[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 2], hm1)

    # Middle state
    ax2 = Axis(fig[1, 3],
        title="Collision (t=T/2)",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm2 = heatmap!(ax2, x, y, rho_tmid[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 4], hm2)

    # Final state
    ax3 = Axis(fig[1, 5],
        title="Final State (t=T)",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm3 = heatmap!(ax3, x, y, rho_final[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 6], hm3)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_soliton_binary(result; save_path=joinpath(FIGURES_DIR, "fig_soliton_binary.png"))

    ψ = result[:ψ]
    x = result[:x]
    y = result[:y]
    z = result[:z]
    Nx = result[:Nx]
    Tmax = result[:Tmax]
    rho_final = abs2.(ψ)
    rho_max = maximum(rho_final)
    iz_center = div(Nx, 2)

    xxx = [x[i] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    yyy = [y[j] for i in 1:Nx, j in 1:Nx, k in 1:Nx]
    zzz = [z[k] for i in 1:Nx, j in 1:Nx, k in 1:Nx]

    rho_t0 = abs2.(analytic_soliton_binary_3d(xxx, yyy, zzz, 0.0))
    rho_tmid = abs2.(analytic_soliton_binary_3d(xxx, yyy, zzz, Tmax/2))

    fig = Figure(size=(1200, 400), fontsize=11)

    # Initial state
    ax1 = Axis(fig[1, 1],
        title="Binary Initial State",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm1 = heatmap!(ax1, x, y, rho_t0[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 2], hm1)

    # Middle state
    ax2 = Axis(fig[1, 3],
        title="Binary at t=T/2",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm2 = heatmap!(ax2, x, y, rho_tmid[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 4], hm2)

    # Final state
    ax3 = Axis(fig[1, 5],
        title="Binary Final State",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)
    hm3 = heatmap!(ax3, x, y, rho_final[:, :, iz_center],
        colormap=:viridis, colorrange=(0, rho_max/3))
    Colorbar(fig[1, 6], hm3)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_tde_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_tde_convergence.png"))

    Nx_list = result[:Nx_list]
    profiles = result[:profiles]
    mass_ratios = result[:mass_ratios]

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: Density profiles overlay
    ax1 = Axis(fig[1, 1],
        title="TDE Density Profile Convergence",
        xlabel="x",
        ylabel="|ψ|²",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red]
    for (i, p) in enumerate(profiles)
        N = length(p.x)
        iy_center = div(N, 2)
        rho_line = p.rho[iy_center, :]
        lines!(ax1, p.x, rho_line, color=colors[i], linewidth=2, label=p.label)
    end
    axislegend(ax1, position=:rt, framevisible=true)

    # Right panel: Mass conservation vs resolution
    ax2 = Axis(fig[1, 2],
        title="TDE Mass Conservation",
        xlabel="Nx",
        ylabel="Mass Remaining (%)",
        xgridvisible=true,
        ygridvisible=true)

    scatter!(ax2, Nx_list, 100 .* mass_ratios,
        color=:blue, marker=:circle, markersize=10)

    # Reference line at 100%
    hlines!(ax2, 100.0, color=:gray, linestyle=:dash, linewidth=1.5, label="100%")

    axislegend(ax2, position=:lb, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    println("\n" * "="^60)
    println("TDE Convergence Results:")
    println("="^60)
    println("Nx\t\tMass Remaining")
    for i in eachindex(Nx_list)
        println("$(Nx_list[i])\t\t$(@sprintf("%.2f%%", 100*mass_ratios[i]))")
    end
    println("="^60)

    return fig
end

function plot_collision_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_collision_convergence.png"))

    Nx_list = result[:Nx_list]
    profiles = result[:profiles]

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: 2D density maps at different resolutions
    ax1 = Axis(fig[1, 1],
        title="Collision Density at Different Resolutions",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red]
    for (i, p) in enumerate(profiles)
        N = length(p.x)
        iz_center = div(N, 2)
        rho_slice = p.rho[:, iz_center]
        hm = heatmap!(ax1, p.x, p.y, rho_slice,
            colormap=:viridis, alpha=0.3)
    end

    # Right panel: Profile line comparison
    ax2 = Axis(fig[1, 2],
        title="Collision Profile Comparison",
        xlabel="x",
        ylabel="|ψ|²",
        xgridvisible=true,
        ygridvisible=true)

    for (i, p) in enumerate(profiles)
        N = length(p.x)
        iz_center = div(N, 2)
        iy_center = div(N, 2)
        rho_line = p.rho[iy_center, :]
        lines!(ax2, p.x, rho_line, color=colors[i], linewidth=2, label=p.label)
    end
    axislegend(ax2, position=:rt, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_binary_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_binary_convergence.png"))

    Nx_list = result[:Nx_list]
    profiles = result[:profiles]

    fig = Figure(size=(1000, 400), fontsize=11)

    # Left panel: 2D density maps at different resolutions
    ax1 = Axis(fig[1, 1],
        title="Binary Density at Different Resolutions",
        xlabel="x",
        ylabel="y",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red]
    for (i, p) in enumerate(profiles)
        N = length(p.x)
        iz_center = div(N, 2)
        rho_slice = p.rho[:, iz_center]
        heatmap!(ax1, p.x, p.y, rho_slice,
            colormap=:viridis, alpha=0.3)
    end

    # Right panel: Profile line comparison
    ax2 = Axis(fig[1, 2],
        title="Binary Profile Comparison",
        xlabel="x",
        ylabel="|ψ|²",
        xgridvisible=true,
        ygridvisible=true)

    for (i, p) in enumerate(profiles)
        N = length(p.x)
        iz_center = div(N, 2)
        iy_center = div(N, 2)
        rho_line = p.rho[iy_center, :]
        lines!(ax2, p.x, rho_line, color=colors[i], linewidth=2, label=p.label)
    end
    axislegend(ax2, position=:rt, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

# ============================================================================
# Main: Run All Tests and Generate Figures
# ============================================================================

@info "="^70
@info "WaveDM.jl - 3D SPE Validation Tests"
@info "Based on PyUltraLight (Edwards et al. 2018)"
@info "="^70

@info "\n[1/7] Running single soliton test..."
result_soliton = run_single_soliton_test(Xmax=30.0, Tmax=2.0, Nx=96, Nt=96)

@info "\n[2/7] Running spatial convergence test..."
result_spatial_conv = run_spatial_convergence_test(Nx_list=[32, 48, 64, 96, 128], Nt=96)

@info "\n[3/7] Running temporal convergence test..."
result_temporal_conv = run_temporal_convergence_test(Nx=96, Nt_list=[32, 48, 64, 96, 128])

@info "\n[4/7] Running soliton collision convergence test..."
result_collision_conv = run_collision_convergence_test(Nx_list=[48, 64, 96, 128])

@info "\n[5/7] Running soliton binary convergence test..."
result_binary_conv = run_binary_convergence_test(Nx_list=[48, 64, 96, 128])

@info "\n[6/7] Running TDE convergence test..."
result_tde_conv = run_tde_convergence_test(Nx_list=[48, 64, 96, 128])

@info "\n[7/7] Generating convergence figures..."
plot_spatial_convergence(result_spatial_conv)
plot_temporal_convergence(result_temporal_conv)
plot_collision_convergence(result_collision_conv)
plot_binary_convergence(result_binary_conv)
plot_tde_convergence(result_tde_conv)

@info "\n" * "="^70
@info "All tests completed!"
@info "Figures saved to: $FIGURES_DIR"
@info "="^70
