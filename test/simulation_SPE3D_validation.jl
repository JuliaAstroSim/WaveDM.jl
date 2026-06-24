#=
3D Schrödinger-Poisson Equation Validation Tests

Usage:
include("test/simulation_SPE3D_validation.jl")
include("E:/JuliaAstroSim/WaveDM.jl/test/simulation_SPE3D_validation.jl")
=#

using WaveDM
using LinearAlgebra
using FFTW
using ProgressMeter
using Statistics
using Printf
using CairoMakie
using FileIO
using Unitful

astro()
const C = Constant(uAstro)

const η₀ = sqrt(C.μ_0 / C.ε_0)
const ħ = C.h/2/π
const Ωₘ₀ = 0.31

const OUTPUT_DIR = "E:/JuliaAstroSim/WaveDM.jl/test/output/SPE3D"
const FIGURES_DIR = joinpath(OUTPUT_DIR, "figures")

mkpath.([OUTPUT_DIR, FIGURES_DIR])

# ============================================================================
# FDM Physical Parameters (same as integration_test.jl)
# ============================================================================
const mₐ = 0.2 * 1e-22 * 1.783e-36 * u"kg"
const aₛ = +8.29e-60 * u"fm"

const length_astro = uconvert(u"kpc", (8 * π * ħ^2 / (3 * mₐ^2 * C.H^2 * Ωₘ₀))^0.25)
const time_astro = uconvert(u"Gyr", (3 * C.H^2 * Ωₘ₀ / (8 * π))^-0.5)
const mass_astro = uconvert(u"Msun", (3 * C.H^2 * Ωₘ₀ / (8 * π))^0.25 * ħ^1.5 / (mₐ^1.5 * C.G))
const velocity_astro = uconvert(u"km/s", length_astro / time_astro)

@info "FDM Physical Units:"
@info "  length_astro = $length_astro"
@info "  time_astro = $time_astro"
@info "  mass_astro = $mass_astro"
@info "  velocity_astro = $velocity_astro"

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
# Test Functions (Physical Units)
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
        Realtime=true,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
        plot_virial = true,
        mₐ=mₐ, aₛ=aₛ, Ωₘ₀=Ωₘ₀,
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
        :Nx => Nx, :Nt => Nt,
        :Xmax => Xmax, :Tmax => Tmax,
        :A => A, :v => v, :x0 => x0,
        :psi_ana => psi_ana,
        :error_l2 => compute_error_l2(ψ_final, psi_ana),
        :error_linf => compute_error_linf(ψ_final, psi_ana),
        :rel_error_l2 => compute_relative_error_l2(ψ_final, psi_ana),
        :dfProp => dfProp
    )
end

function run_soliton_collision_test(; Xmax=60.0, Tmax=100.0, Nx=64, Nt=64, verbose=true)

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
        Realtime=true,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
        plot_virial = true,
        mₐ=mₐ, aₛ=aₛ, Ωₘ₀=Ωₘ₀,
    )

    return Dict(
        :ψ => ψ_final,
        :x => collect(LinRange(-Xmax, Xmax, Nx)),
        :y => collect(LinRange(-Xmax, Xmax, Nx)),
        :z => collect(LinRange(-Xmax, Xmax, Nx)),
        :Nx => Nx, :Nt => Nt,
        :Xmax => Xmax, :Tmax => Tmax,
        :dfProp => dfProp
    )
end

function analytic_soliton_tde_3d(x, y, z, t; A=sqrt(2.0), x0=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0])
    return analytic_soliton_3d(x, y, z, t; A=A, v=v, x0=x0)
end

function run_tde_test(; Xmax=50.0, Tmax=100.0, Nx=64, Nt=64,
    A=sqrt(2.0),
    soliton_x0=[5.0, 5.0, 0.0],
    soliton_v=[0.0, 0.0, 0.0],
    perturber_mass=2.0,
    perturber_x0=[8.0, 8.0, 0.0],
    perturber_v=[-0.3, -0.3, 0.0],
    softening=1.0,
    verbose=true)

    verbose && @info "Running TDE test: M_perturber=$perturber_mass"

    function IC_with_perturber(x, y, z)
        soliton = analytic_soliton_3d(x, y, z, 0.0; A=A, v=soliton_v, x0=soliton_x0)
        return soliton
    end

    V_tde = let perturber_x0=perturber_x0, perturber_mass=perturber_mass, softening=softening
        (x, y, z, ψ) -> begin
            r_perturber = sqrt.((x .- perturber_x0[1]).^2 + (y .- perturber_x0[2]).^2 + (z .- perturber_x0[3]).^2)
            G = 1.0
            M = perturber_mass
            potential = -G * M ./ sqrt.(r_perturber.^2 .+ softening^2)
            return potential
        end
    end

    ψ_final, fig, figMOND, chi2RC, dfProp, dfAcc, averaged_ψ2 = SPE3D_waveDM(;
        Xmax=Xmax, Ymax=Xmax, Zmax=Xmax, Tmax=Tmax,
        Nx=Nx, Ny=Nx, Nz=Nx, Nt=Nt,
        IC=IC_with_perturber,
        V=V_tde,
        baryon_mode=:ignored,
        absorb_coeff=0.0,
        outputdir=joinpath(OUTPUT_DIR, "tde"),
        title="TDE",
        gpu=false,
        Realtime=true,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
        plot_virial = true,
        mₐ=mₐ, aₛ=aₛ, Ωₘ₀=Ωₘ₀,
    )

    return Dict(
        :ψ => ψ_final,
        :x => collect(LinRange(-Xmax, Xmax, Nx)),
        :y => collect(LinRange(-Xmax, Xmax, Nx)),
        :z => collect(LinRange(-Xmax, Xmax, Nx)),
        :Nx => Nx, :Nt => Nt,
        :Xmax => Xmax, :Tmax => Tmax,
        :perturber_mass => perturber_mass,
        :dfProp => dfProp
    )
end

function run_soliton_binary_test(; Xmax=50.0, Tmax=300.0, Nx=64, Nt=64,
    A=sqrt(2.0), separation=20.0, orbital_velocity=0.2,
    verbose=true)

    verbose && @info "Running soliton binary test: Nx=$Nx, Nt=$Nt"

    IC = (x, y, z) -> analytic_soliton_binary_3d(x, y, z, 0.0;
        A=A, separation=separation, orbital_velocity=orbital_velocity)

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
        Realtime=true,
        plotOptical=false,
        plotWaveDM=false,
        FDM_mass_ratio=1.0,
        KDK_flag=true,
        best_fit_halo_mass=false,
        save_IC=false,
        save_phi=false,
        plot_virial = true,
        mₐ=mₐ, aₛ=aₛ, Ωₘ₀=Ωₘ₀,
    )

    return Dict(
        :ψ => ψ_final,
        :x => collect(LinRange(-Xmax, Xmax, Nx)),
        :y => collect(LinRange(-Xmax, Xmax, Nx)),
        :z => collect(LinRange(-Xmax, Xmax, Nx)),
        :Nx => Nx, :Nt => Nt,
        :Xmax => Xmax, :Tmax => Tmax,
        :dfProp => dfProp
    )
end

function run_tde_convergence_test(; Xmax=50.0, Tmax=100.0, Nx_list=[48, 64, 96], Nt_list=nothing, verbose=true)
    if isnothing(Nt_list)
        Nt_list = Nx_list
    end

    verbose && @info "Running TDE convergence test"

    profiles = []
    results = []

    for (i, Nx) in enumerate(Nx_list)
        verbose && @info "  Nx = $Nx"
        Nt = Nt_list[min(i, length(Nt_list))]
        result = run_tde_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        push!(results, result)

        x = result[:x]
        iz_center = div(Nx, 2)
        rho_slice = abs2.(result[:ψ])[:, :, iz_center]
        push!(profiles, (x=x, y=x, rho=rho_slice, label="Nx=$Nx"))
    end

    return Dict(
        :Nx_list => Nx_list,
        :profiles => profiles,
        :results => results
    )
end

function run_collision_convergence_test(; Xmax=60.0, Tmax=100.0, Nx_list=[48, 64, 96], Nt_list=nothing, verbose=true)
    if isnothing(Nt_list)
        Nt_list = Nx_list
    end

    verbose && @info "Running collision convergence test"

    profiles = []
    results = []

    for (i, Nx) in enumerate(Nx_list)
        verbose && @info "  Nx = $Nx"
        Nt = Nt_list[min(i, length(Nt_list))]
        result = run_soliton_collision_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        x = result[:x]
        N = length(x)
        iy_center = div(N, 2)
        iz_center = div(N, 2)
        rho_line = abs2.(result[:ψ])[iy_center, iz_center, :]
        push!(profiles, (x=x, rho=rho_line, label="Nx=$Nx"))
        push!(results, result)
    end

    return Dict(
        :Nx_list => Nx_list,
        :profiles => profiles,
        :results => results
    )
end

function run_binary_convergence_test(; Xmax=50.0, Tmax=300.0, Nx_list=[48, 64, 96], Nt_list=nothing, verbose=true)
    if isnothing(Nt_list)
        Nt_list = Nx_list
    end

    verbose && @info "Running binary convergence test"

    profiles = []
    results = []

    for (i, Nx) in enumerate(Nx_list)
        verbose && @info "  Nx = $Nx"
        Nt = Nt_list[min(i, length(Nt_list))]
        result = run_soliton_binary_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        x = result[:x]
        N = length(x)
        iy_center = div(N, 2)
        iz_center = div(N, 2)
        rho_line = abs2.(result[:ψ])[iy_center, iz_center, :]
        push!(profiles, (x=x, rho=rho_line, label="Nx=$Nx"))
        push!(results, result)
    end

    return Dict(
        :Nx_list => Nx_list,
        :profiles => profiles,
        :results => results
    )
end

# ============================================================================
# Temporal Convergence Tests for Collision, Binary, and TDE
# ============================================================================

function run_collision_temporal_convergence_test(; Xmax=60.0, Tmax=100.0, Nx=64, Nt_list=[48, 64, 96], verbose=true)
    verbose && @info "Running collision temporal convergence test"

    profiles = []
    results = []

    for Nt in Nt_list
        verbose && @info "  Nt = $Nt"
        result = run_soliton_collision_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        x = result[:x]
        iy_center = div(Nx, 2)
        iz_center = div(Nx, 2)
        rho_line = abs2.(result[:ψ])[iy_center, iz_center, :]
        push!(profiles, (x=x, rho=rho_line, label="Nt=$Nt"))
        push!(results, result)
    end

    return Dict(
        :Nt_list => Nt_list,
        :profiles => profiles,
        :results => results
    )
end

function run_binary_temporal_convergence_test(; Xmax=50.0, Tmax=300.0, Nx=64, Nt_list=[48, 64, 96], verbose=true)
    verbose && @info "Running binary temporal convergence test"

    profiles = []
    results = []

    for Nt in Nt_list
        verbose && @info "  Nt = $Nt"
        result = run_soliton_binary_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)

        x = result[:x]
        iy_center = div(Nx, 2)
        iz_center = div(Nx, 2)
        rho_line = abs2.(result[:ψ])[iy_center, iz_center, :]
        push!(profiles, (x=x, rho=rho_line, label="Nt=$Nt"))
        push!(results, result)
    end

    return Dict(
        :Nt_list => Nt_list,
        :profiles => profiles,
        :results => results
    )
end

function run_tde_temporal_convergence_test(; Xmax=50.0, Tmax=100.0, Nx=64, Nt_list=[48, 64, 96], verbose=true)
    verbose && @info "Running TDE temporal convergence test"

    results = []

    for Nt in Nt_list
        verbose && @info "  Nt = $Nt"
        result = run_tde_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)
        push!(results, result)
    end

    return Dict(
        :Nt_list => Nt_list,
        :results => results
    )
end

function run_spatial_convergence_test(; Xmax=30.0, Tmax=2.0, Nx_list=[32, 64], Nt=96, verbose=true)

    verbose && @info "Running spatial convergence test"

    errors_l2 = Float64[]
    errors_linf = Float64[]
    rel_errors_l2 = Float64[]
    profiles = []
    results = []

    for Nx in Nx_list
        verbose && @info "  Nx = $Nx"
        result = run_single_soliton_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)
        push!(errors_l2, result[:error_l2])
        push!(errors_linf, result[:error_linf])
        push!(rel_errors_l2, result[:rel_error_l2])
        push!(profiles, (x=result[:x], rho=abs2.(result[:ψ]), label="Nx=$Nx"))
        push!(results, result)
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
        :profiles => profiles,
        :results => results
    )
end

function run_temporal_convergence_test(; Xmax=30.0, Tmax=2.0, Nx=96, Nt_list=[32, 64], verbose=true)

    verbose && @info "Running temporal convergence test"

    errors_l2 = Float64[]
    errors_linf = Float64[]
    rel_errors_l2 = Float64[]
    profiles = []
    results = []

    for Nt in Nt_list
        verbose && @info "  Nt = $Nt"
        result = run_single_soliton_test(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, verbose=false)
        push!(errors_l2, result[:error_l2])
        push!(errors_linf, result[:error_linf])
        push!(rel_errors_l2, result[:rel_error_l2])
        push!(profiles, (x=result[:x], rho=abs2.(result[:ψ]), label="Nt=$Nt"))
        push!(results, result)
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
        :profiles => profiles,
        :results => results  # Include full results for dfProp
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
    profiles = result[:profiles]

    fig = Figure(size=(700, 500), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    ax = Axis(fig[1, 1],
        title="Single Soliton - Spatial Resolution Convergence",
        xlabel="r",
        ylabel="|ψ|²",
        xgridvisible=false,
        ygridvisible=false)

    for (i, p) in enumerate(profiles)
        Nx = length(p.x)
        iz_center = div(Nx, 2)
        iy_center = div(Nx, 2)
        rho_line = p.rho[iy_center, iz_center, :]
        r_line = abs.(p.x)
        label = "$(Nx_list[i])³"
        lines!(ax, r_line, rho_line, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax, position=:rt, framevisible=true)
    xlims!(ax, 0, 5)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_temporal_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_temporal_convergence.png"))

    Nt_list = result[:Nt_list]
    profiles = result[:profiles]

    fig = Figure(size=(700, 500), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    ax = Axis(fig[1, 1],
        title="Single Soliton - Temporal Resolution Convergence",
        xlabel="r",
        ylabel="|ψ|²",
        xgridvisible=false,
        ygridvisible=false)

    for (i, p) in enumerate(profiles)
        Nx = length(p.x)
        iz_center = div(Nx, 2)
        iy_center = div(Nx, 2)
        rho_line = p.rho[iy_center, iz_center, :]
        r_line = abs.(p.x)
        label = "Nt=$(Nt_list[i])"
        lines!(ax, r_line, rho_line, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax, position=:rt, framevisible=true)
    xlims!(ax, 0, 5)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

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

function plot_tde_convergence(spatial_result, temporal_result=nothing; save_path=joinpath(FIGURES_DIR, "fig_tde_convergence.png"))

    spatial_results = spatial_result[:results]
    spatial_Nx_list = spatial_result[:Nx_list]

    fig = Figure(size=(1000, 800), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    # Top panel: Spatial resolution - Energy conservation over time
    ax_top = Axis(fig[1, 1],
        title="TDE - Spatial Resolution Convergence",
        xlabel="Time",
        ylabel=L"\Delta E_{total} / E_{total}",
        xgridvisible=false,
        ygridvisible=false)

    for (i, res) in enumerate(spatial_results)
        if !haskey(res, :dfProp) || isnothing(res[:dfProp])
            continue
        end
        dfProp = res[:dfProp]
        if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
            continue
        end

        t = dfProp.t
        PE = dfProp.PE_abs
        KE = dfProp.KE
        QE = dfProp.QE
        E_total = KE .+ QE .- PE

        if length(E_total) > 0 && E_total[1] != 0
            dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
            label = "$(spatial_Nx_list[i])³"
            lines!(ax_top, t, dE_E, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
    end
    hlines!(ax_top, [0], color=:black, linestyle=:dash)
    axislegend(ax_top, position=:rt, framevisible=true)

    # Bottom panel: Temporal resolution - Energy conservation over time
    if !isnothing(temporal_result)
        temporal_results = temporal_result[:results]
        temporal_Nt_list = temporal_result[:Nt_list]

        ax_bottom = Axis(fig[2, 1],
            title="TDE - Temporal Resolution Convergence",
            xlabel="Time",
            ylabel=L"\Delta E_{total} / E_{total}",
            xgridvisible=false,
            ygridvisible=false)

        for (i, res) in enumerate(temporal_results)
            if !haskey(res, :dfProp) || isnothing(res[:dfProp])
                continue
            end
            dfProp = res[:dfProp]
            if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
                continue
            end

            t = dfProp.t
            PE = dfProp.PE_abs
            KE = dfProp.KE
            QE = dfProp.QE
            E_total = KE .+ QE .- PE

            if length(E_total) > 0 && E_total[1] != 0
                dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
                label = "Nt=$(temporal_Nt_list[i])"
                lines!(ax_bottom, t, dE_E, color=colors[min(i, length(colors))],
                       linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
            end
        end
        hlines!(ax_bottom, [0], color=:black, linestyle=:dash)
        axislegend(ax_bottom, position=:rt, framevisible=true)
    else
        # Fallback: Just show final energy vs resolution
        ax_bottom = Axis(fig[2, 1],
            title="TDE - Final Energy vs Resolution",
            xlabel="Grid Size Nx",
            ylabel="Total Energy",
            xgridvisible=false,
            ygridvisible=false)

        total_energies = Float64[]
        for res in spatial_results
            if !haskey(res, :dfProp) || isnothing(res[:dfProp])
                push!(total_energies, NaN)
                continue
            end
            dfProp = res[:dfProp]
            if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
                push!(total_energies, NaN)
                continue
            end
            PE = dfProp.PE_abs
            KE = dfProp.KE
            QE = dfProp.QE
            E_total = KE .+ QE .- PE
            push!(total_energies, E_total[1])
        end
        scatterlines!(ax_bottom, spatial_Nx_list, total_energies,
            color=:blue, linewidth=2, marker=:circle, markersize=10)
    end

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_collision_convergence_energy(spatial_result, temporal_result=nothing; save_path=joinpath(FIGURES_DIR, "fig_collision_convergence_energy.png"))

    spatial_results = spatial_result[:results]
    spatial_Nx_list = spatial_result[:Nx_list]

    fig = Figure(size=(1000, 800), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    ax_top = Axis(fig[1, 1],
        title="Collision - Spatial Resolution Convergence",
        xlabel="Time",
        ylabel=L"\Delta E_{total} / E_{total}",
        xgridvisible=false,
        ygridvisible=false)

    for (i, res) in enumerate(spatial_results)
        if !haskey(res, :dfProp) || isnothing(res[:dfProp])
            continue
        end
        dfProp = res[:dfProp]
        if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
            continue
        end

        t = dfProp.t
        PE = dfProp.PE_abs
        KE = dfProp.KE
        QE = dfProp.QE
        E_total = KE .+ QE .- PE

        if length(E_total) > 0 && E_total[1] != 0
            dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
            label = "$(spatial_Nx_list[i])³"
            lines!(ax_top, t, dE_E, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
    end
    hlines!(ax_top, [0], color=:black, linestyle=:dash)
    axislegend(ax_top, position=:rt, framevisible=true)

    if !isnothing(temporal_result)
        temporal_results = temporal_result[:results]
        temporal_Nt_list = temporal_result[:Nt_list]

        ax_bottom = Axis(fig[2, 1],
            title="Collision - Temporal Resolution Convergence",
            xlabel="Time",
            ylabel=L"\Delta E_{total} / E_{total}",
            xgridvisible=false,
            ygridvisible=false)

        for (i, res) in enumerate(temporal_results)
            if !haskey(res, :dfProp) || isnothing(res[:dfProp])
                continue
            end
            dfProp = res[:dfProp]
            if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
                continue
            end

            t = dfProp.t
            PE = dfProp.PE_abs
            KE = dfProp.KE
            QE = dfProp.QE
            E_total = KE .+ QE .- PE

            if length(E_total) > 0 && E_total[1] != 0
                dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
                label = "Nt=$(temporal_Nt_list[i])"
                lines!(ax_bottom, t, dE_E, color=colors[min(i, length(colors))],
                       linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
            end
        end
        hlines!(ax_bottom, [0], color=:black, linestyle=:dash)
        axislegend(ax_bottom, position=:rt, framevisible=true)
    else
        ax_bottom = Axis(fig[2, 1],
            title="Collision - Final Energy vs Resolution",
            xlabel="Grid Size Nx",
            ylabel="Total Energy",
            xgridvisible=false,
            ygridvisible=false)

        total_energies = Float64[]
        for res in spatial_results
            if !haskey(res, :dfProp) || isnothing(res[:dfProp])
                push!(total_energies, NaN)
                continue
            end
            dfProp = res[:dfProp]
            if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
                push!(total_energies, NaN)
                continue
            end
            PE = dfProp.PE_abs
            KE = dfProp.KE
            QE = dfProp.QE
            E_total = KE .+ QE .- PE
            push!(total_energies, E_total[1])
        end
        scatterlines!(ax_bottom, spatial_Nx_list, total_energies,
            color=:blue, linewidth=2, marker=:circle, markersize=10)
    end

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_binary_convergence_energy(spatial_result, temporal_result=nothing; save_path=joinpath(FIGURES_DIR, "fig_binary_convergence_energy.png"))

    spatial_results = spatial_result[:results]
    spatial_Nx_list = spatial_result[:Nx_list]

    fig = Figure(size=(1000, 800), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    ax_top = Axis(fig[1, 1],
        title="Binary - Spatial Resolution Convergence",
        xlabel="Time",
        ylabel=L"\Delta E_{total} / E_{total}",
        xgridvisible=false,
        ygridvisible=false)

    for (i, res) in enumerate(spatial_results)
        if !haskey(res, :dfProp) || isnothing(res[:dfProp])
            continue
        end
        dfProp = res[:dfProp]
        if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
            continue
        end

        t = dfProp.t
        PE = dfProp.PE_abs
        KE = dfProp.KE
        QE = dfProp.QE
        E_total = KE .+ QE .- PE

        if length(E_total) > 0 && E_total[1] != 0
            dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
            label = "$(spatial_Nx_list[i])³"
            lines!(ax_top, t, dE_E, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
    end
    hlines!(ax_top, [0], color=:black, linestyle=:dash)
    axislegend(ax_top, position=:rt, framevisible=true)

    if !isnothing(temporal_result)
        temporal_results = temporal_result[:results]
        temporal_Nt_list = temporal_result[:Nt_list]

        ax_bottom = Axis(fig[2, 1],
            title="Binary - Temporal Resolution Convergence",
            xlabel="Time",
            ylabel=L"\Delta E_{total} / E_{total}",
            xgridvisible=false,
            ygridvisible=false)

        for (i, res) in enumerate(temporal_results)
            if !haskey(res, :dfProp) || isnothing(res[:dfProp])
                continue
            end
            dfProp = res[:dfProp]
            if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
                continue
            end

            t = dfProp.t
            PE = dfProp.PE_abs
            KE = dfProp.KE
            QE = dfProp.QE
            E_total = KE .+ QE .- PE

            if length(E_total) > 0 && E_total[1] != 0
                dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
                label = "Nt=$(temporal_Nt_list[i])"
                lines!(ax_bottom, t, dE_E, color=colors[min(i, length(colors))],
                       linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
            end
        end
        hlines!(ax_bottom, [0], color=:black, linestyle=:dash)
        axislegend(ax_bottom, position=:rt, framevisible=true)
    else
        ax_bottom = Axis(fig[2, 1],
            title="Binary - Final Energy vs Resolution",
            xlabel="Grid Size Nx",
            ylabel="Total Energy",
            xgridvisible=false,
            ygridvisible=false)

        total_energies = Float64[]
        for res in spatial_results
            if !haskey(res, :dfProp) || isnothing(res[:dfProp])
                push!(total_energies, NaN)
                continue
            end
            dfProp = res[:dfProp]
            if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
                push!(total_energies, NaN)
                continue
            end
            PE = dfProp.PE_abs
            KE = dfProp.KE
            QE = dfProp.QE
            E_total = KE .+ QE .- PE
            push!(total_energies, E_total[1])
        end
        scatterlines!(ax_bottom, spatial_Nx_list, total_energies,
            color=:blue, linewidth=2, marker=:circle, markersize=10)
    end

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_collision_convergence(spatial_result, temporal_result=nothing; save_path=joinpath(FIGURES_DIR, "fig_collision_convergence.png"))

    spatial_profiles = spatial_result[:profiles]
    spatial_Nx_list = spatial_result[:Nx_list]

    fig = Figure(size=(700, 700), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    ax_top = Axis(fig[1, 1],
        title="Collision - Spatial Resolution Convergence",
        xlabel="x",
        ylabel="Density",
        xgridvisible=false,
        ygridvisible=false)

    for (i, p) in enumerate(spatial_profiles)
        label = "$(spatial_Nx_list[i])³"
        lines!(ax_top, p.x, p.rho, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax_top, position=:rt, framevisible=true)
    xlims!(ax_top, -10, 10)

    if !isnothing(temporal_result)
        temporal_profiles = temporal_result[:profiles]
        temporal_Nt_list = temporal_result[:Nt_list]

        ax_bottom = Axis(fig[2, 1],
            title="Collision - Temporal Resolution Convergence",
            xlabel="x",
            ylabel="Density",
            xgridvisible=false,
            ygridvisible=false)

        for (i, p) in enumerate(temporal_profiles)
            label = "Nt=$(temporal_Nt_list[i])"
            lines!(ax_bottom, p.x, p.rho, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
        axislegend(ax_bottom, position=:rt, framevisible=true)
    else
        ax_bottom = Axis(fig[2, 1],
            title="Collision - Density Profile Detail",
            xlabel="x",
            ylabel="Density",
            xgridvisible=false,
            ygridvisible=false)

        for (i, p) in enumerate(spatial_profiles)
            label = "$(spatial_Nx_list[i])³"
            lines!(ax_bottom, p.x, p.rho, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
        axislegend(ax_bottom, position=:rt, framevisible=true)
        xlims!(ax_bottom, -10, 10)
    end

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_binary_convergence(spatial_result, temporal_result=nothing; save_path=joinpath(FIGURES_DIR, "fig_binary_convergence.png"))

    spatial_profiles = spatial_result[:profiles]
    spatial_Nx_list = spatial_result[:Nx_list]

    fig = Figure(size=(700, 700), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :solid, :dash, :dashdot]

    ax_top = Axis(fig[1, 1],
        title="Binary System - Spatial Resolution Convergence",
        xlabel="x",
        ylabel="Density",
        xgridvisible=false,
        ygridvisible=false)

    for (i, p) in enumerate(spatial_profiles)
        label = "$(spatial_Nx_list[i])³"
        lines!(ax_top, p.x, p.rho, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax_top, position=:rt, framevisible=true)
    xlims!(ax_top, -15, 15)

    if !isnothing(temporal_result)
        temporal_profiles = temporal_result[:profiles]
        temporal_Nt_list = temporal_result[:Nt_list]

        ax_bottom = Axis(fig[2, 1],
            title="Binary System - Temporal Resolution Convergence",
            xlabel="x",
            ylabel="Density",
            xgridvisible=false,
            ygridvisible=false)

        for (i, p) in enumerate(temporal_profiles)
            label = "Nt=$(temporal_Nt_list[i])"
            lines!(ax_bottom, p.x, p.rho, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
        axislegend(ax_bottom, position=:rt, framevisible=true)
    else
        ax_bottom = Axis(fig[2, 1],
            title="Binary System - Density Profile Detail",
            xlabel="x",
            ylabel="Density",
            xgridvisible=false,
            ygridvisible=false)

        for (i, p) in enumerate(spatial_profiles)
            label = "$(spatial_Nx_list[i])³"
            lines!(ax_bottom, p.x, p.rho, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
        axislegend(ax_bottom, position=:rt, framevisible=true)
        xlims!(ax_bottom, -15, 15)
    end

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

# ============================================================================
# Temporal Convergence Plotting Functions
# ============================================================================

function plot_collision_temporal_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_collision_temporal_convergence.png"))
    Nt_list = result[:Nt_list]
    profiles = result[:profiles]

    fig = Figure(size=(700, 500), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :dash, :dashdot]

    ax = Axis(fig[1, 1],
        title="Collision - Temporal Resolution Convergence",
        xlabel="x",
        ylabel="Density",
        xgridvisible=false,
        ygridvisible=false)

    for (i, p) in enumerate(profiles)
        label = "Nt=$(Nt_list[i])"
        lines!(ax, p.x, p.rho, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax, position=:rt, framevisible=true)
    xlims!(ax, -10, 10)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_binary_temporal_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_binary_temporal_convergence.png"))
    Nt_list = result[:Nt_list]
    profiles = result[:profiles]

    fig = Figure(size=(700, 500), fontsize=12)

    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :solid, :dash, :dashdot]

    ax = Axis(fig[1, 1],
        title="Binary System - Temporal Resolution Convergence",
        xlabel="x",
        ylabel="Density",
        xgridvisible=false,
        ygridvisible=false)

    for (i, p) in enumerate(profiles)
        label = "Nt=$(Nt_list[i])"
        lines!(ax, p.x, p.rho, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax, position=:rt, framevisible=true)
    xlims!(ax, -15, 15)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_tde_temporal_convergence(result; save_path=joinpath(FIGURES_DIR, "fig_tde_temporal_convergence.png"))
    Nt_list = result[:Nt_list]
    results = result[:results]

    fig = Figure(size=(1000, 400), fontsize=12)

    ax1 = Axis(fig[1, 1],
        title="TDE - Energy Convergence",
        xlabel="Time",
        ylabel="Total Energy",
        xgridvisible=false,
        ygridvisible=false)

    ax2 = Axis(fig[1, 2],
        title="TDE - Energy Conservation Error",
        xlabel="Time",
        ylabel=L"\Delta E / E_0",
        xgridvisible=false,
        ygridvisible=false)

    colors = [:blue, :orange, :green, :red, :purple]

    for (i, res) in enumerate(results)
        if !haskey(res, :dfProp) || isnothing(res[:dfProp])
            continue
        end
        dfProp = res[:dfProp]
        if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
            continue
        end

        t = dfProp.t
        PE = dfProp.PE_abs
        KE = dfProp.KE
        QE = dfProp.QE
        E_total = KE .+ QE .- PE

        label = "Nt=$(Nt_list[i])"
        lines!(ax1, t, E_total, color=colors[min(i, length(colors))], linewidth=2, label=label)

        if length(E_total) > 0 && E_total[1] != 0
            dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
            lines!(ax2, t, dE_E, color=colors[min(i, length(colors))], linewidth=2, label=label)
        end
    end

    axislegend(ax1, position=:rt, framevisible=true)
    axislegend(ax2, position=:rt, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

# ============================================================================
# Virial Energy Terms Plotting
# Based on the reference figure showing energy evolution
# ============================================================================

function plot_virial_energy_terms(dfProp; save_path=joinpath(FIGURES_DIR, "fig_virial_energy.png"), title="Virial Energy Terms")
    if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
        @warn "dfProp missing virial energy columns (PE_abs, KE, QE). Skipping virial energy plot."
        return nothing
    end

    t = dfProp.t
    PE = dfProp.PE_abs  # Potential energy (absolute value)
    KE = dfProp.KE      # Kinetic energy
    QE = dfProp.QE      # Quantum energy

    # Energy components as shown in the reference figure
    E_K_plus_Q = KE .+ QE           # Kinetic + Quantum (orange line in ref)
    E_GP = -PE                      # Gravitational potential (green line in ref, negative because PE is stored as absolute)
    E_total = KE .+ QE .- PE        # Total energy (blue line in ref)

    fig = Figure(size=(800, 600), fontsize=14)

    ax = Axis(fig[1, 1],
        title=title,
        xlabel="Time (code units)",
        ylabel="Energy (code units)",
        xgridvisible=true,
        ygridvisible=true)

    # Plot energy components matching the reference figure style
    # E_K + E_Q (orange line)
    lines!(ax, t, E_K_plus_Q, color=:orange, linewidth=2.5, label=L"E_K + E_Q")

    # Total energy (blue line)
    lines!(ax, t, E_total, color=:blue, linewidth=2.5, label="Total energy")

    # E_GP self-interaction (green line)
    lines!(ax, t, E_GP, color=:green, linewidth=2.5, label=L"E_{GP} (self-interaction)")

    # Add text annotations similar to reference figure
    text!(ax, t[1] + 0.05*(t[end]-t[1]), maximum(E_K_plus_Q)*0.9, text=L"E_K + E_Q", fontsize=12, color=:orange)
    text!(ax, t[1] + 0.05*(t[end]-t[1]), E_total[1]*0.9, text="Total energy", fontsize=12, color=:blue)
    text!(ax, t[1] + 0.05*(t[end]-t[1]), minimum(E_GP)*0.9, text=L"E_{GP} (self-interaction)", fontsize=12, color=:green)

    axislegend(ax, position=:rt, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_virial_energy_convergence(results::Vector{Dict}; save_path=joinpath(FIGURES_DIR, "fig_virial_energy_convergence.png"), title="Virial Energy Convergence")
    fig = Figure(size=(1000, 400), fontsize=12)

    # Left panel: Total energy comparison
    ax1 = Axis(fig[1, 1],
        title="Total Energy Convergence",
        xlabel="Time (code units)",
        ylabel="Total Energy",
        xgridvisible=true,
        ygridvisible=true)

    # Right panel: Energy components for finest resolution
    ax2 = Axis(fig[1, 2],
        title="Energy Components (Finest Resolution)",
        xlabel="Time (code units)",
        ylabel="Energy (code units)",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red, :purple]

    for (i, result) in enumerate(results)
        if !haskey(result, :dfProp) || isnothing(result[:dfProp])
            continue
        end
        dfProp = result[:dfProp]
        if !hasproperty(dfProp, :PE_abs) || !hasproperty(dfProp, :KE) || !hasproperty(dfProp, :QE)
            continue
        end

        t = dfProp.t
        PE = dfProp.PE_abs
        KE = dfProp.KE
        QE = dfProp.QE
        E_total = KE .+ QE .- PE

        label = haskey(result, :Nx) ? "Nx=$(result[:Nx])" : "Run $i"
        lines!(ax1, t, E_total, color=colors[min(i, length(colors))], linewidth=2, label=label)

        # Only plot components for the last (finest resolution) result
        if i == length(results)
            lines!(ax2, t, KE .+ QE, color=:orange, linewidth=2.5, label=L"E_K + E_Q")
            lines!(ax2, t, E_total, color=:blue, linewidth=2.5, label="Total energy")
            lines!(ax2, t, -PE, color=:green, linewidth=2.5, label=L"E_{GP}")
        end
    end

    axislegend(ax1, position=:rt, framevisible=true)
    axislegend(ax2, position=:rt, framevisible=true)

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

# ============================================================================
# Momentum Conservation Plot
# Based on Prop.jl showing momentum components over time
# ============================================================================

function plot_momentum_conservation(dfProp; save_path=joinpath(FIGURES_DIR, "fig_momentum_conservation.png"), title="Momentum Conservation")
    if !hasproperty(dfProp, :MomentumX) || !hasproperty(dfProp, :MomentumY) || !hasproperty(dfProp, :MomentumZ)
        @warn "dfProp missing momentum columns (MomentumX, MomentumY, MomentumZ). Skipping momentum plot."
        return nothing
    end

    t = dfProp.t
    px = dfProp.MomentumX
    py = dfProp.MomentumY
    pz = dfProp.MomentumZ

    fig = Figure(size=(1600, 600), fontsize=22)

    ax = Axis(fig[1, 1],
        title=title,
        xlabel="t [Gyr]",
        ylabel="p [Msun*kpc/Gyr]",
        xgridvisible=true,
        ygridvisible=true)

    line_px = lines!(ax, t, px, color=:red, linewidth=2, label=L"p_x")
    line_py = lines!(ax, t, py, color=:green, linewidth=2, label=L"p_y")
    line_pz = lines!(ax, t, pz, color=:blue, linewidth=2, label=L"p_z")

    hlines!(ax, [0]; color=:black, linestyle=:dash)
    xlims!(ax, extrema(t))

    Legend(fig[1, 2],
        [line_px, line_py, line_pz],
        [L"p_x", L"p_y", L"p_z"];
        tellheight=false,
        margin=(10, 10, 10, 10),
        backgroundcolor=(:white, 0.5),
    )

    save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_momentum_conservation_convergence(results::Vector{Dict}; save_path=joinpath(FIGURES_DIR, "fig_momentum_conservation_convergence.png"), title="Momentum Conservation Convergence")
    fig = Figure(size=(1000, 400), fontsize=12)

    # Left panel: Momentum magnitude comparison
    ax1 = Axis(fig[1, 1],
        title="Momentum Magnitude Convergence",
        xlabel="Time (code units)",
        ylabel="|p| [Msun*kpc/Gyr]",
        xgridvisible=true,
        ygridvisible=true)

    # Right panel: Momentum components for finest resolution
    ax2 = Axis(fig[1, 2],
        title="Momentum Components (Finest Resolution)",
        xlabel="Time (code units)",
        ylabel="p [Msun*kpc/Gyr]",
        xgridvisible=true,
        ygridvisible=true)

    colors = [:blue, :orange, :green, :red, :purple]

    for (i, result) in enumerate(results)
        if !haskey(result, :dfProp) || isnothing(result[:dfProp])
            continue
        end
        dfProp = result[:dfProp]
        if !hasproperty(dfProp, :MomentumX) || !hasproperty(dfProp, :MomentumY) || !hasproperty(dfProp, :MomentumZ)
            continue
        end

        t = dfProp.t
        px = dfProp.MomentumX
        py = dfProp.MomentumY
        pz = dfProp.MomentumZ
        p_mag = sqrt.(px.^2 .+ py.^2 .+ pz.^2)

        label = haskey(result, :Nx) ? "$(result[:Nx])³" : "Run $i"
        lines!(ax1, t, p_mag, color=colors[min(i, length(colors))], linewidth=2, label=label)

        # Only plot components for the last (finest resolution) result
        if i == length(results)
            lines!(ax2, t, px, color=:red, linewidth=2, label=L"p_x")
            lines!(ax2, t, py, color=:green, linewidth=2, label=L"p_y")
            lines!(ax2, t, pz, color=:blue, linewidth=2, label=L"p_z")
            hlines!(ax2, [0]; color=:black, linestyle=:dash)
        end
    end

    axislegend(ax1, position=:rt, framevisible=true)
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

@info "\n[1/9] Running single soliton spatial convergence test..."
result_spatial_conv = run_spatial_convergence_test(Nx_list=[32, 64], Nt=96)

@info "\n[2/9] Running single soliton temporal convergence test..."
result_temporal_conv = run_temporal_convergence_test(Nx=96, Nt_list=[32, 64])

@info "\n[3/9] Running soliton collision spatial convergence test..."
result_collision_conv = run_collision_convergence_test(Nx_list=[32, 64])

@info "\n[4/9] Running soliton collision temporal convergence test..."
result_collision_temporal = run_collision_temporal_convergence_test(Nx=64, Nt_list=[32, 64])

@info "\n[5/9] Running soliton binary spatial convergence test..."
result_binary_conv = run_binary_convergence_test(Nx_list=[32, 64])

@info "\n[6/9] Running soliton binary temporal convergence test..."
result_binary_temporal = run_binary_temporal_convergence_test(Nx=64, Nt_list=[32, 64])

@info "\n[7/9] Running TDE spatial convergence test..."
result_tde_conv = run_tde_convergence_test(Nx_list=[32, 64])

@info "\n[8/9] Running TDE temporal convergence test..."
result_tde_temporal = run_tde_temporal_convergence_test(Nx=64, Nt_list=[32, 64])

@info "\n[9/9] Generating all figures..."

# Spatial convergence plots
plot_spatial_convergence(result_spatial_conv)
plot_collision_convergence(result_collision_conv, result_collision_temporal)
plot_collision_convergence_energy(result_collision_conv, result_collision_temporal)
plot_binary_convergence(result_binary_conv, result_binary_temporal)
plot_binary_convergence_energy(result_binary_conv, result_binary_temporal)
plot_tde_convergence(result_tde_conv, result_tde_temporal)

# Temporal convergence plots
plot_temporal_convergence(result_temporal_conv)
plot_collision_temporal_convergence(result_collision_temporal)
plot_binary_temporal_convergence(result_binary_temporal)
plot_tde_temporal_convergence(result_tde_temporal)

# Momentum conservation (use finest resolution)
if !isempty(result_spatial_conv[:results])
    finest_result = result_spatial_conv[:results][end]
    if haskey(finest_result, :dfProp) && !isnothing(finest_result[:dfProp])
        plot_momentum_conservation(finest_result[:dfProp],
            save_path=joinpath(FIGURES_DIR, "fig_momentum_conservation_single_soliton.png"),
            title="Momentum Conservation - Single Soliton")
    end
end

@info "\n" * "="^70
@info "All tests completed!"
@info "Figures saved to: $FIGURES_DIR"
@info "="^70
