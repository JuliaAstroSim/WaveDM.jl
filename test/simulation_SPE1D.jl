using LinearAlgebra
using FFTW
using ProgressMeter
using Statistics
using Printf
using CairoMakie

function analytic_harmonic_oscillator(x, t)
    return exp.(-x.^2 ./ 2) .* exp.(-im .* t ./ 2)
end

function analytic_soliton(x, t; A=sqrt(2.0), v=0.0, x0=0.0)
    xi = x .- x0 .- v .* t
    B = A / 2
    omega = A^2 / 2 - v^2 / 2
    k = v
    return A .* sech.(B .* xi) .* exp.(im .* (k .* x .- omega .* t))
end

function analytic_soliton_collision(x, t; A1=sqrt(2.0), v1=0.5, A2=sqrt(2.0), v2=-0.5, x1=10.0, x2=-10.0)
    sol1 = analytic_soliton(x, t; A=A1, v=v1, x0=x1)
    sol2 = analytic_soliton(x, t; A=A2, v=v2, x0=x2)
    return sol1 .+ sol2
end

function compute_error_l2(psi_num, psi_ana)
    return sqrt(mean(abs.(psi_num .- psi_ana).^2))
end

function compute_error_linf(psi_num, psi_ana)
    return maximum(abs.(psi_num .- psi_ana))
end

function compute_energy(ψ, x, Δ, V_func)
    Nx = length(x)
    dx = Δ
    
    # Kinetic energy: E_K = ∫ |∇ψ|²/2 dx
    k = fftfreq(Nx, 2π/dx)
    ψ_hat = fft(ψ)
    grad_ψ_hat = im .* k .* ψ_hat
    grad_ψ = ifft(grad_ψ_hat)
    E_K = sum(abs2.(grad_ψ)) * dx / 2
    
    # Quantum pressure: E_Q = ∫ |∇|ψ||²/2 dx
    rho = abs2.(ψ)
    sqrt_rho = sqrt.(rho)
    grad_sqrt_rho_hat = im .* k .* fft(sqrt_rho)
    grad_sqrt_rho = ifft(grad_sqrt_rho_hat)
    E_Q = sum(abs2.(grad_sqrt_rho)) * dx / 2
    
    # External potential energy: E_ext = ∫ V|ψ|² dx
    E_ext = sum(V_func.(x, ψ) .* rho) * dx
    
    # Self-interaction (Gross-Pitaevskii): E_GP = -∫ |ψ|⁴/2 dx
    E_GP = -sum(rho.^2) * dx / 2
    
    E_total = real(E_K + E_Q + E_ext + E_GP)
    
    return Dict(:E_K=>real(E_K), :E_Q=>real(E_Q), :E_K_plus_E_Q=>real(E_K+E_Q), 
                :E_ext=>real(E_ext), :E_GP=>real(E_GP), :E_total=>E_total)
end

function setup_fft_operators_1d(Xmax, Nx, dt)
    kx = collect(LinRange(-Nx/4/Xmax, Nx/4/Xmax-1/2/Xmax, Nx))
    Laplacian = [(2π*im*kx[i])^2 for i in 1:Nx]
    linear_phase = fftshift(exp.(im * Laplacian * dt / 2))
    return linear_phase
end

function setup_absorption_boundary_1d(Xmax, absorb_coeff, x, dt)
    wx = Xmax/40
    border = exp.(-absorb_coeff .* (2 .- tanh.((x.+Xmax)./wx) .+ tanh.((x.-Xmax)./wx)) .* dt)
    return border
end

function SPE1D(;
    Xmax=20.0,
    Tmax=1.0,
    Nx=256,
    Nt=128,
    IC=x->exp.(-x.^2 ./ 2),
    V=(x, ψ)->0.0,
    absorb_coeff=0.0,
    compute_energy_flag=false,
    analytic_solution=nothing,
    outputdir="./output/SPE1D",
    filename="SPE1D",
    title="1D SPE",
    verbose=true
)
    verbose && @info "Initializing 1D SPE simulation: $title"

    mkpath(outputdir)

    x = collect(LinRange(-Xmax, Xmax, Nx))
    t = collect(LinRange(0, Tmax, Nt))
    dt = Tmax / Nt
    Δ = x[2] - x[1]

    Uout = zeros(ComplexF64, Nx, Nt)
    Uout[:, 1] = IC.(x)

    linear_phase = setup_fft_operators_1d(Xmax, Nx, dt)
    border = setup_absorption_boundary_1d(Xmax, absorb_coeff, x, dt)

    # Energy tracking
    if compute_energy_flag
        E_K = zeros(Nt)
        E_Q = zeros(Nt)
        E_K_plus_E_Q = zeros(Nt)
        E_GP = zeros(Nt)
        E_total = zeros(Nt)
        
        energy_dict = compute_energy(Uout[:, 1], x, Δ, V)
        E_K[1] = energy_dict[:E_K]
        E_Q[1] = energy_dict[:E_Q]
        E_K_plus_E_Q[1] = energy_dict[:E_K_plus_E_Q]
        E_GP[1] = energy_dict[:E_GP]
        E_total[1] = energy_dict[:E_total]
    end

    if verbose
        @showprogress "Propagating..." for i in 2:Nt
            ψ = Uout[:, i-1]
            potential = V.(x, ψ)
            spec = fft(ψ .* exp.(-im * dt * potential))
            Uout[:, i] = border .* ifft(spec .* linear_phase)
            
            if compute_energy_flag
                energy_dict = compute_energy(Uout[:, i], x, Δ, V)
                E_K[i] = energy_dict[:E_K]
                E_Q[i] = energy_dict[:E_Q]
                E_K_plus_E_Q[i] = energy_dict[:E_K_plus_E_Q]
                E_GP[i] = energy_dict[:E_GP]
                E_total[i] = energy_dict[:E_total]
            end
        end
    else
        for i in 2:Nt
            ψ = Uout[:, i-1]
            potential = V.(x, ψ)
            spec = fft(ψ .* exp.(-im * dt * potential))
            Uout[:, i] = border .* ifft(spec .* linear_phase)
            
            if compute_energy_flag
                energy_dict = compute_energy(Uout[:, i], x, Δ, V)
                E_K[i] = energy_dict[:E_K]
                E_Q[i] = energy_dict[:E_Q]
                E_K_plus_E_Q[i] = energy_dict[:E_K_plus_E_Q]
                E_GP[i] = energy_dict[:E_GP]
                E_total[i] = energy_dict[:E_total]
            end
        end
    end

    result = Dict{Symbol, Any}()
    result[:ψ] = Uout
    result[:x] = x
    result[:t] = t
    result[:dt] = dt
    result[:Δ] = Δ
    result[:Nx] = Nx
    result[:Nt] = Nt
    result[:Xmax] = Xmax
    result[:Tmax] = Tmax
    result[:outputdir] = outputdir
    result[:filename] = filename
    result[:title] = title

    if compute_energy_flag
        result[:E_K] = E_K
        result[:E_Q] = E_Q
        result[:E_K_plus_E_Q] = E_K_plus_E_Q
        result[:E_GP] = E_GP
        result[:E_total] = E_total
        result[:delta_E_total] = (E_total .- E_total[1]) ./ abs(E_total[1])
    end

    if analytic_solution !== nothing
        error_l2 = zeros(Nt)
        error_linf = zeros(Nt)
        for i in 1:Nt
            psi_ana = analytic_solution(x, t[i])
            error_l2[i] = compute_error_l2(Uout[:, i], psi_ana)
            error_linf[i] = compute_error_linf(Uout[:, i], psi_ana)
        end
        result[:error_l2] = error_l2
        result[:error_linf] = error_linf
        result[:error_l2_final] = error_l2[end]
        result[:error_linf_final] = error_linf[end]
    end

    return result
end

function test_harmonic_oscillator(; Xmax=20.0, Tmax=2.0, Nx=512, Nt=256, outputdir="./output/SPE1D/harmonic_oscillator", verbose=true)
    verbose && @info "Running Harmonic Oscillator Test"

    IC = x->exp.(-x.^2 ./ 2)
    analytic_sol(x, t) = analytic_harmonic_oscillator(x, t)
    V = (x, ψ)->x.^2 ./ 2

    result = SPE1D(; Xmax, Tmax, Nx, Nt, IC, V, absorb_coeff=0.0,
        analytic_solution=analytic_sol, outputdir, filename="harmonic_oscillator",
        title="Harmonic Oscillator (1D)", verbose)

    verbose && @info "Final L2 Error: $(result[:error_l2_final])"
    verbose && @info "Final L∞ Error: $(result[:error_linf_final])"

    return result
end

function test_single_soliton(; Xmax=30.0, Tmax=3.0, Nx=512, Nt=256, A=sqrt(2.0), v=0.0, x0=0.0, outputdir="./output/SPE1D/single_soliton", verbose=true)
    verbose && @info "Running Single Soliton Test (A=$A, v=$v, x0=$x0)"

    IC = x->analytic_soliton(x, 0.0; A, v, x0)
    analytic_sol(x, t) = analytic_soliton(x, t; A, v, x0)
    V = (x, ψ)->-abs.(ψ).^2

    result = SPE1D(; Xmax, Tmax, Nx, Nt, IC, V, absorb_coeff=0.0,
        analytic_solution=analytic_sol, outputdir, filename="single_soliton_A$(A)_v$(v)",
        title="Single Soliton (1D, A=$A, v=$v)", verbose)

    verbose && @info "Final L2 Error: $(result[:error_l2_final])"
    verbose && @info "Final L∞ Error: $(result[:error_linf_final])"

    return result
end

function test_soliton_collision(; Xmax=60.0, Tmax=10.0, Nx=512, Nt=512, A1=sqrt(2.0), v1=0.5, A2=sqrt(2.0), v2=-0.5, x1=10.0, x2=-10.0, outputdir="./output/SPE1D/soliton_collision", verbose=true)
    verbose && @info "Running Soliton Collision Test (v1=$v1, v2=$v2)"

    IC = x->analytic_soliton_collision(x, 0.0; A1, v1, A2, v2, x1, x2)
    analytic_sol(x, t) = analytic_soliton_collision(x, t; A1, v1, A2, v2, x1, x2)
    V = (x, ψ)->-abs.(ψ).^2

    result = SPE1D(; Xmax, Tmax, Nx, Nt, IC, V, absorb_coeff=0.0,
        analytic_solution=analytic_sol, outputdir, filename="soliton_collision",
        title="Soliton Collision (1D)", verbose)

    verbose && @info "Final L2 Error: $(result[:error_l2_final])"
    verbose && @info "Final L∞ Error: $(result[:error_linf_final])"

    return result
end

function convergence_spatial(; Xmax=20.0, Tmax=1.0, Nx_list=[64, 128, 256, 512, 1024], Nt=512, test_case="harmonic_oscillator", outputdir="./output/SPE1D/convergence", verbose=true)
    verbose && @info "Running Spatial Convergence Test for $test_case"

    error_l2 = Float64[]
    error_linf = Float64[]

    for Nx in Nx_list
        if test_case == "harmonic_oscillator"
            result = test_harmonic_oscillator(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, outputdir=outputdir, verbose=false)
        elseif test_case == "single_soliton"
            result = test_single_soliton(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, outputdir=outputdir, verbose=false)
        else
            error("Unknown test case: $test_case")
        end

        push!(error_l2, result[:error_l2_final])
        push!(error_linf, result[:error_linf_final])
    end

    dx_list = 2 * Xmax ./ Nx_list
    order_l2 = Float64[]
    order_linf = Float64[]

    for i in 2:length(Nx_list)
        push!(order_l2, log(error_l2[i] / error_l2[i-1]) / log(dx_list[i] / dx_list[i-1]))
        push!(order_linf, log(error_linf[i] / error_linf[i-1]) / log(dx_list[i] / dx_list[i-1]))
    end

    if verbose
        @info "Spatial Convergence Results for $test_case:"
        @info "Nx\t\tL2 Error\t\tL∞ Error\t\tOrder L2\t\tOrder L∞"
        for i in 1:length(Nx_list)
            order_l2_str = i > 1 ? @sprintf("%.2f", order_l2[i-1]) : "-"
            order_linf_str = i > 1 ? @sprintf("%.2f", order_linf[i-1]) : "-"
            @info "$(Nx_list[i])\t\t$(error_l2[i])\t\t$(error_linf[i])\t\t$order_l2_str\t\t$order_linf_str"
        end
    end

    return Dict(:Nx_list=>Nx_list, :error_l2=>error_l2, :error_linf=>error_linf, :order_l2=>order_l2, :order_linf=>order_linf, :test_case=>test_case)
end

function convergence_temporal(; Xmax=20.0, Tmax=1.0, Nx=512, Nt_list=[32, 64, 128, 256, 512], test_case="harmonic_oscillator", outputdir="./output/SPE1D/convergence", verbose=true)
    verbose && @info "Running Temporal Convergence Test for $test_case"

    error_l2 = Float64[]
    error_linf = Float64[]

    for Nt in Nt_list
        if test_case == "harmonic_oscillator"
            result = test_harmonic_oscillator(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, outputdir=outputdir, verbose=false)
        elseif test_case == "single_soliton"
            result = test_single_soliton(Xmax=Xmax, Tmax=Tmax, Nx=Nx, Nt=Nt, outputdir=outputdir, verbose=false)
        else
            error("Unknown test case: $test_case")
        end

        push!(error_l2, result[:error_l2_final])
        push!(error_linf, result[:error_linf_final])
    end

    dt_list = Tmax ./ Nt_list
    order_l2 = Float64[]
    order_linf = Float64[]

    for i in 2:length(Nt_list)
        push!(order_l2, log(error_l2[i] / error_l2[i-1]) / log(dt_list[i] / dt_list[i-1]))
        push!(order_linf, log(error_linf[i] / error_linf[i-1]) / log(dt_list[i] / dt_list[i-1]))
    end

    if verbose
        @info "Temporal Convergence Results for $test_case:"
        @info "Nt\t\tL2 Error\t\tL∞ Error\t\tOrder L2\t\tOrder L∞"
        for i in 1:length(Nt_list)
            order_l2_str = i > 1 ? @sprintf("%.2f", order_l2[i-1]) : "-"
            order_linf_str = i > 1 ? @sprintf("%.2f", order_linf[i-1]) : "-"
            @info "$(Nt_list[i])\t\t$(error_l2[i])\t\t$(error_linf[i])\t\t$order_l2_str\t\t$order_linf_str"
        end
    end

    return Dict(:Nt_list=>Nt_list, :error_l2=>error_l2, :error_linf=>error_linf, :order_l2=>order_l2, :order_linf=>order_linf, :test_case=>test_case)
end

# ==================== Plotting Functions ====================

function plot_energy_conservation(results_dict::Dict; outputdir="./output/SPE1D/plots", filename="energy_conservation.png")
    mkpath(outputdir)
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
        xlabel="Time (code units)", 
        ylabel="ΔE_total / |E_total|",
        title="Energy Conservation vs Resolution")
    
    colors = [:blue, :orange, :green, :red, :purple]
    
    for (i, (label, result)) in enumerate(results_dict)
        if haskey(result, :delta_E_total)
            t = result[:t]
            delta_E = result[:delta_E_total]
            lines!(ax, t, delta_E .* 1e5, label=label, color=colors[mod1(i, length(colors))], linewidth=2)
        end
    end
    
    axislegend(ax, position=:rt)
    save(joinpath(outputdir, filename), fig)
    @info "Saved energy conservation plot to $(joinpath(outputdir, filename))"
    return fig
end

function plot_energy_components(result::Dict; outputdir="./output/SPE1D/plots", filename="energy_components.png")
    mkpath(outputdir)
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
        xlabel="Time (code units)", 
        ylabel="Energy (code units)",
        title="Energy Components Evolution")
    
    t = result[:t]
    
    lines!(ax, t, result[:E_K_plus_E_Q], label="E_K + E_Q", color=:orange, linewidth=2)
    lines!(ax, t, result[:E_GP], label="E_GP (self-interaction)", color=:green, linewidth=2)
    lines!(ax, t, result[:E_total], label="Total energy", color=:blue, linewidth=2)
    
    axislegend(ax, position=:lt)
    save(joinpath(outputdir, filename), fig)
    @info "Saved energy components plot to $(joinpath(outputdir, filename))"
    return fig
end

function plot_density_scaling(results_dict::Dict; outputdir="./output/SPE1D/plots", filename="density_scaling.png", time_idx=nothing)
    mkpath(outputdir)
    
    # Determine layout based on number of results
    n_results = length(results_dict)
    if n_results <= 3
        fig = Figure(size=(800, 500))
        ax = Axis(fig[1, 1], xlabel="Position", ylabel="Density", title="Density Profile vs Resolution")
    else
        fig = Figure(size=(800, 800))
        ax1 = Axis(fig[1, 1], xlabel="Position", ylabel="Density", title="Spatial Resolution Scaling")
        ax2 = Axis(fig[2, 1], xlabel="Position", ylabel="Density", title="Temporal Resolution Scaling")
    end
    
    colors = [:blue, :orange, :green, :red, :purple, :brown]
    linestyles = [nothing, nothing, nothing, :dash, :dashdot, :dot]
    
    # Plot spatial resolution scaling
    for (i, (label, result)) in enumerate(results_dict)
        x = result[:x]
        psi = result[:ψ]
        t = result[:t]
        
        # Use middle time point if not specified
        idx = time_idx === nothing ? div(size(psi, 2), 2) : time_idx
        
        rho = abs2.(psi[:, idx])
        
        if n_results <= 3
            lines!(ax, x, rho, label=label, color=colors[mod1(i, length(colors))], 
                   linestyle=linestyles[mod1(i, length(linestyles))], linewidth=2)
        else
            if i <= 3
                lines!(ax1, x, rho, label=label, color=colors[mod1(i, length(colors))], linewidth=2)
            else
                lines!(ax2, x, rho, label=label, color=colors[mod1(i, length(colors))], linewidth=2)
            end
        end
    end
    
    if n_results <= 3
        axislegend(ax, position=:rt)
    else
        axislegend(ax1, position=:rt)
        axislegend(ax2, position=:rt)
    end
    
    save(joinpath(outputdir, filename), fig)
    @info "Saved density scaling plot to $(joinpath(outputdir, filename))"
    return fig
end

function run_scaling_convergence_tests(; verbose=true)
    verbose && @info "=" ^ 60
    verbose && @info "Running Scaling Convergence Tests with Energy Tracking"
    verbose && @info "=" ^ 60
    
    outputdir = "./output/SPE1D/scaling"
    mkpath(outputdir)
    
    # Test 1: Soliton collision with different spatial resolutions (like Figure 6 & 11 top)
    verbose && @info "\n--- Soliton Collision: Spatial Resolution Scaling ---"
    Nx_list = [256, 512, 1024]
    Nt = 1024
    Tmax = 6.0
    
    spatial_results = Dict{String, Dict}()
    for Nx in Nx_list
        label = "$(Nx)"
        verbose && @info "Running with Nx=$Nx, Nt=$Nt"
        result = SPE1D(
            Xmax=60.0, Tmax=Tmax, Nx=Nx, Nt=Nt,
            IC=x->analytic_soliton_collision(x, 0.0; A1=sqrt(2.0), v1=0.3, A2=sqrt(2.0), v2=-0.3, x1=15.0, x2=-15.0),
            V=(x,ψ)->-abs.(ψ).^2,
            compute_energy_flag=true,
            absorb_coeff=0.0,
            outputdir=outputdir,
            filename="soliton_collision_Nx$(Nx)",
            title="Soliton Collision Nx=$Nx",
            verbose=false
        )
        spatial_results[label] = result
    end
    
    # Plot energy conservation (Figure 6 style)
    plot_energy_conservation(spatial_results, outputdir=outputdir, filename="fig6_energy_conservation_spatial.png")
    
    # Plot density scaling (Figure 11 top style)
    plot_density_scaling(spatial_results, outputdir=outputdir, filename="fig11_density_spatial.png")
    
    # Test 2: Soliton collision with different temporal resolutions (like Figure 11 bottom)
    verbose && @info "\n--- Soliton Collision: Temporal Resolution Scaling ---"
    Nx = 512
    Nt_list = [256, 512, 1024, 2048]
    Tmax = 6.0
    
    temporal_results = Dict{String, Dict}()
    for Nt in Nt_list
        label = "$(Nt)"
        verbose && @info "Running with Nx=$Nx, Nt=$Nt"
        result = SPE1D(
            Xmax=60.0, Tmax=Tmax, Nx=Nx, Nt=Nt,
            IC=x->analytic_soliton_collision(x, 0.0; A1=sqrt(2.0), v1=0.3, A2=sqrt(2.0), v2=-0.3, x1=15.0, x2=-15.0),
            V=(x,ψ)->-abs.(ψ).^2,
            compute_energy_flag=true,
            absorb_coeff=0.0,
            outputdir=outputdir,
            filename="soliton_collision_Nt$(Nt)",
            title="Soliton Collision Nt=$Nt",
            verbose=false
        )
        temporal_results[label] = result
    end
    
    # Plot energy conservation for temporal scaling
    plot_energy_conservation(temporal_results, outputdir=outputdir, filename="fig6_energy_conservation_temporal.png")
    
    # Plot density scaling for temporal (Figure 11 bottom style)
    plot_density_scaling(temporal_results, outputdir=outputdir, filename="fig11_density_temporal.png")
    
    # Test 3: Energy components evolution (Figure 5 style)
    verbose && @info "\n--- Energy Components Evolution ---"
    result_energy = SPE1D(
        Xmax=60.0, Tmax=6.0, Nx=512, Nt=1024,
        IC=x->analytic_soliton_collision(x, 0.0; A1=sqrt(2.0), v1=0.3, A2=sqrt(2.0), v2=-0.3, x1=15.0, x2=-15.0),
        V=(x,ψ)->-abs.(ψ).^2,
        compute_energy_flag=true,
        absorb_coeff=0.0,
        outputdir=outputdir,
        filename="soliton_collision_energy",
        title="Soliton Collision Energy",
        verbose=false
    )
    
    plot_energy_components(result_energy, outputdir=outputdir, filename="fig5_energy_components.png")
    
    verbose && @info "\n" * "=" ^ 60
    verbose && @info "Scaling Convergence Tests Complete"
    verbose && @info "Plots saved to: $outputdir"
    verbose && @info "=" ^ 60
    
    return Dict(:spatial=>spatial_results, :temporal=>temporal_results, :energy=>result_energy)
end

function run_all_1D_SPE_tests(; verbose=true)
    verbose && @info "=" ^ 60
    verbose && @info "Running All 1D SPE Validation Tests"
    verbose && @info "=" ^ 60

    results = Dict{String, Any}()

    verbose && @info "\n--- Test 1: Harmonic Oscillator ---"
    results["harmonic_oscillator"] = test_harmonic_oscillator(Nx=512, Nt=256, verbose=verbose)

    verbose && @info "\n--- Test 2: Single Soliton ---"
    results["single_soliton"] = test_single_soliton(Nx=512, Nt=256, verbose=verbose)

    verbose && @info "\n--- Test 3: Soliton Collision ---"
    results["soliton_collision"] = test_soliton_collision(Nx=512, Nt=512, verbose=verbose)

    verbose && @info "\n--- Convergence Analysis: Spatial (Harmonic Oscillator) ---"
    results["convergence_spatial"] = convergence_spatial(Nx_list=[64, 128, 256, 512], Nt=512, test_case="harmonic_oscillator", verbose=verbose)

    verbose && @info "\n--- Convergence Analysis: Temporal (Harmonic Oscillator) ---"
    results["convergence_temporal"] = convergence_temporal(Xmax=20.0, Tmax=1.0, Nx=512, Nt_list=[32, 64, 128, 256, 512], test_case="harmonic_oscillator", verbose=verbose)

    verbose && @info "\n" * "=" ^ 60
    verbose && @info "Test Summary"
    verbose && @info "=" ^ 60

    for (name, result) in results
        if name in ["harmonic_oscillator", "single_soliton", "soliton_collision"]
            if haskey(result, :error_l2_final)
                verbose && @info "$name: L2 Error = $(result[:error_l2_final]), L∞ Error = $(result[:error_linf_final])"
            end
        end
    end

    return results
end

if basename(pwd()) == "test" || isfile(joinpath(@__DIR__, "simulation_SPE1D.jl"))
    @info "Running 1D SPE tests..."
    
    # Run basic tests
    # results = run_all_1D_SPE_tests()
    
    # Run scaling convergence tests with plotting
    @info "\nRunning scaling convergence tests with energy tracking and plotting..."
    scaling_results = run_scaling_convergence_tests()
    
else
    @info "To run tests, include this file from the test directory or call run_all_1D_SPE_tests()"
end
