#=
Energy Convergence Tests for WaveDM Prop Output
通用化的能量收敛测试画图函数

Usage:
includet("test/energy_convergence.jl")

Nx_files = ["path/to/Nx64.csv", "path/to/Nx128.csv", "path/to/Nx256.csv"]
Nt_files = ["path/to/Nt58.csv", "path/to/Nt116.csv", ...]
plot_energy_convergence(Nx_files, Nt_files; outputdir="output/path")
=#

using WaveDM
using DataFrames
using CSV
using CairoMakie
using FileIO
using Interpolations

function compute_energy_conservation(df)
    required_cols = [:PE_abs, :KE, :QE]
    all(hasproperty(df, c) for c in required_cols) || return nothing, nothing, nothing, nothing

    t = df.t
    PE = df.PE_abs
    KE = df.KE
    QE = df.QE
    E_total = KE .+ QE .- PE

    return t, E_total, KE, QE
end

function compute_mass_conservation(df)
    hasproperty(df, :TotalMass) || return nothing, nothing
    
    t = df.t
    M = df.TotalMass
    
    return t, M
end

function interpolate_to_common_grid(t_target, t_source, E_source)
    itp = interpolate((t_source,), E_source, Gridded(Linear()))
    return itp.(t_target)
end

function parse_Nx_Nt(filepath::AbstractString)
    m1 = match(r"Nx\((\d+)\)", filepath)
    m2 = match(r"Nt\((\d+)\)", filepath)
    return isnothing(m1) ? nothing : parse(Int, m1[1]),
        isnothing(m2) ? nothing : parse(Int, m2[1])
end

function plot_energy_convergence(Nx_files, Nt_files;
    fig_size::Tuple{Int, Int}=(800, 1000),
    save_prefix::AbstractString="fig_energy_convergence",
    xlabel::AbstractString="T [Gyr]",
    ylabel::AbstractString=L"\Delta E_{total} / E_{total}",
    spatial_label_fn::Function=(Nx, Nt) -> "$(Nx)³ (Nt=$Nt)",
    temporal_label_fn::Union{Function, Nothing}=nothing,
    temporal_nt_ref::Union{Int, Nothing}=nothing
)
    figures_dir = joinpath(outputdir, "figures")
    mkpath(figures_dir)

    spatial_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nx_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(spatial_data, (df, Nx, Nt))
    end
    sort!(spatial_data, by=x->x[2])
    spatial_ref_E = spatial_data[end][1] |> compute_energy_conservation |> x->x[2]

    temporal_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nt_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(temporal_data, (df, Nx, Nt))
    end
    sort!(temporal_data, by=x->x[3])
    temporal_ref_E = temporal_data[end][1] |> compute_energy_conservation |> x->x[2]

    isempty(spatial_data) && isempty(temporal_data) && (@warn "No valid data loaded") && return nothing

    @info "Plotting energy convergence"
    fig = Figure(size=fig_size, fontsize=14)
    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :dash, :dashdot, :dot, :solid]

    ax_top = Axis(fig[1, 1],
        title="Spatial Resolution Convergence (ref: $(spatial_data[end][2])³)",
        xlabel=xlabel,
        ylabel=L"E - E_{ref}",
        xgridvisible=false,
        ygridvisible=false)

    t_ref = spatial_data[end][1].t
    spatial_ref_E_interp = spatial_ref_E

    for (i, (df, Nx, Nt)) in enumerate(spatial_data)
        t, E_total, _, _ = compute_energy_conservation(df)
        isnothing(t) && continue
        length(E_total) < 2 && continue

        E_interp = interpolate_to_common_grid(t_ref, t, E_total)
        dE = E_interp .- spatial_ref_E_interp
        
        if Nx == spatial_data[end][2]
            label = "$(Nx)³ (ref)"
            lines!(ax_top, t_ref, dE, color=:black, linewidth=2, linestyle=:dash, label=label)
        else
            label = spatial_label_fn(Nx, Nt)
            lines!(ax_top, t_ref, dE, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
    end
    hlines!(ax_top, [0], color=:black, linestyle=:dash)
    axislegend(ax_top, position=:rt, framevisible=true)

    ax_bottom = Axis(fig[2, 1],
        title="Temporal Resolution Convergence",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    if !isempty(temporal_data)
        Nx_fixed = temporal_data[1][2]
        ax_bottom.title = "Temporal Resolution Convergence (Nx=$Nx_fixed, ref: $(temporal_data[end][3]))"
    end

    t_ref_temporal = temporal_data[end][1].t
    temporal_ref_E_interp = temporal_ref_E

    default_temporal_label_fn = let nt_ref=temporal_nt_ref
        if isnothing(nt_ref)
            Nt -> "Nt=$Nt"
        else
            function(nt::Int)
                ratio = nt_ref / nt
                if abs(ratio - 1) < 1e-6
                    "Δt ×1.0"
                else
                    "Δt ×$(round(ratio, digits=1))"
                end
            end
        end
    end
    actual_temporal_label_fn = isnothing(temporal_label_fn) ? default_temporal_label_fn : temporal_label_fn

    for (i, (df, Nx, Nt)) in enumerate(temporal_data)
        t, E_total, _, _ = compute_energy_conservation(df)
        isnothing(t) && continue
        length(E_total) < 2 && continue

        E_interp = interpolate_to_common_grid(t_ref_temporal, t, E_total)
        dE = E_interp .- temporal_ref_E_interp
        
        if Nt == temporal_data[end][3]
            label = "Nt=$Nt (ref)"
            lines!(ax_bottom, t_ref_temporal, dE, color=:black, linewidth=2, linestyle=:dash, label=label)
        else
            label = actual_temporal_label_fn(Nt)
            lines!(ax_bottom, t_ref_temporal, dE, color=colors[min(i, length(colors))],
                   linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
        end
    end
    hlines!(ax_bottom, [0], color=:black, linestyle=:dash)
    axislegend(ax_bottom, position=:rt, framevisible=true)

    save_path = joinpath(figures_dir, "$save_prefix.png")
    Makie.save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_mass_convergence(Nx_files, Nt_files;
    fig_size::Tuple{Int, Int}=(800, 1000),
    save_prefix::AbstractString="fig_mass_convergence",
    xlabel::AbstractString="T [Gyr]",
    ylabel::AbstractString=L"\Delta M / M_{initial}",
    spatial_label_fn::Function=(Nx, Nt) -> "$(Nx)³ (Nt=$Nt)",
    temporal_label_fn::Union{Function, Nothing}=nothing,
    temporal_nt_ref::Union{Int, Nothing}=nothing
)
    figures_dir = joinpath(outputdir, "figures")
    mkpath(figures_dir)

    spatial_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nx_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(spatial_data, (df, Nx, Nt))
    end
    sort!(spatial_data, by=x->x[2])

    temporal_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nt_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(temporal_data, (df, Nx, Nt))
    end
    sort!(temporal_data, by=x->x[3])

    isempty(spatial_data) && isempty(temporal_data) && (@warn "No valid data loaded") && return nothing

    @info "Plotting mass convergence"
    fig = Figure(size=fig_size, fontsize=14)
    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :dash, :dashdot, :dot, :solid]

    ax_top = Axis(fig[1, 1],
        title="Spatial Resolution - Mass Conservation",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    for (i, (df, Nx, Nt)) in enumerate(spatial_data)
        t, M = compute_mass_conservation(df)
        isnothing(t) && continue
        length(M) < 2 && continue
        M[1] == 0 && continue

        dM_M = (M .- M[1]) ./ M[1]
        label = spatial_label_fn(Nx, Nt)
        lines!(ax_top, t, dM_M, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    hlines!(ax_top, [0], color=:black, linestyle=:dash)
    axislegend(ax_top, position=:lt, framevisible=true)

    ax_bottom = Axis(fig[2, 1],
        title="Temporal Resolution - Mass Conservation",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    if !isempty(temporal_data)
        Nx_fixed = temporal_data[1][2]
        ax_bottom.title = "Temporal Resolution - Mass Conservation (Nx=$Nx_fixed)"
    end

    default_temporal_label_fn = let nt_ref=temporal_nt_ref
        if isnothing(nt_ref)
            Nt -> "Nt=$Nt"
        else
            function(nt::Int)
                ratio = nt_ref / nt
                if abs(ratio - 1) < 1e-6
                    "Δt ×1.0"
                else
                    "Δt ×$(round(ratio, digits=1))"
                end
            end
        end
    end
    actual_temporal_label_fn = isnothing(temporal_label_fn) ? default_temporal_label_fn : temporal_label_fn

    for (i, (df, Nx, Nt)) in enumerate(temporal_data)
        t, M = compute_mass_conservation(df)
        isnothing(t) && continue
        length(M) < 2 && continue
        M[1] == 0 && continue

        dM_M = (M .- M[1]) ./ M[1]
        label = actual_temporal_label_fn(Nt)
        lines!(ax_bottom, t, dM_M, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    hlines!(ax_bottom, [0], color=:black, linestyle=:dash)
    axislegend(ax_bottom, position=:lt, framevisible=true)

    save_path = joinpath(figures_dir, "$save_prefix.png")
    Makie.save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_energy_absolute(Nx_files, Nt_files;
    fig_size::Tuple{Int, Int}=(800, 1000),
    save_prefix::AbstractString="fig_energy_absolute",
    xlabel::AbstractString="T [Gyr]",
    ylabel::AbstractString="E_total",
    spatial_label_fn::Function=(Nx, Nt) -> "$(Nx)³ (Nt=$Nt)",
    temporal_label_fn::Union{Function, Nothing}=nothing,
    temporal_nt_ref::Union{Int, Nothing}=nothing
)
    figures_dir = joinpath(outputdir, "figures")
    mkpath(figures_dir)

    spatial_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nx_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(spatial_data, (df, Nx, Nt))
    end
    sort!(spatial_data, by=x->x[2])

    temporal_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nt_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(temporal_data, (df, Nx, Nt))
    end
    sort!(temporal_data, by=x->x[3])

    isempty(spatial_data) && isempty(temporal_data) && (@warn "No valid data loaded") && return nothing

    @info "Plotting absolute energy"
    fig = Figure(size=fig_size, fontsize=14)
    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :dash, :dashdot, :dot, :solid]

    ax_top = Axis(fig[1, 1],
        title="Spatial Resolution - Absolute Energy",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    for (i, (df, Nx, Nt)) in enumerate(spatial_data)
        t, E_total, _, _ = compute_energy_conservation(df)
        isnothing(t) && continue
        length(E_total) < 2 && continue

        label = spatial_label_fn(Nx, Nt)
        lines!(ax_top, t, E_total, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax_top, position=:rt, framevisible=true)

    ax_bottom = Axis(fig[2, 1],
        title="Temporal Resolution - Absolute Energy",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    if !isempty(temporal_data)
        Nx_fixed = temporal_data[1][2]
        ax_bottom.title = "Temporal Resolution - Absolute Energy (Nx=$Nx_fixed)"
    end

    default_temporal_label_fn = let nt_ref=temporal_nt_ref
        if isnothing(nt_ref)
            Nt -> "Nt=$Nt"
        else
            function(nt::Int)
                ratio = nt_ref / nt
                if abs(ratio - 1) < 1e-6
                    "Δt ×1.0"
                else
                    "Δt ×$(round(ratio, digits=1))"
                end
            end
        end
    end
    actual_temporal_label_fn = isnothing(temporal_label_fn) ? default_temporal_label_fn : temporal_label_fn

    for (i, (df, Nx, Nt)) in enumerate(temporal_data)
        t, E_total, _, _ = compute_energy_conservation(df)
        isnothing(t) && continue
        length(E_total) < 2 && continue

        label = actual_temporal_label_fn(Nt)
        lines!(ax_bottom, t, E_total, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    axislegend(ax_bottom, position=:rt, framevisible=true)

    save_path = joinpath(figures_dir, "$save_prefix.png")
    Makie.save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

function plot_energy_relative(Nx_files, Nt_files;
    fig_size::Tuple{Int, Int}=(800, 1000),
    save_prefix::AbstractString="fig_energy_relative",
    xlabel::AbstractString="T [Gyr]",
    ylabel::AbstractString=L"\Delta E / E_{initial}",
    spatial_label_fn::Function=(Nx, Nt) -> "$(Nx)³ (Nt=$Nt)",
    temporal_label_fn::Union{Function, Nothing}=nothing,
    temporal_nt_ref::Union{Int, Nothing}=nothing
)
    figures_dir = joinpath(outputdir, "figures")
    mkpath(figures_dir)

    spatial_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nx_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(spatial_data, (df, Nx, Nt))
    end
    sort!(spatial_data, by=x->x[2])

    temporal_data = Tuple{DataFrame, Int, Int}[]
    for filepath in Nt_files
        Nx, Nt = parse_Nx_Nt(filepath)
        (isnothing(Nx) || isnothing(Nt)) && continue
        df = CSV.read(filepath, DataFrame)
        push!(temporal_data, (df, Nx, Nt))
    end
    sort!(temporal_data, by=x->x[3])

    isempty(spatial_data) && isempty(temporal_data) && (@warn "No valid data loaded") && return nothing

    @info "Plotting relative energy"
    fig = Figure(size=fig_size, fontsize=14)
    colors = [:blue, :orange, :green, :red, :purple]
    linestyles = [:solid, :dash, :dashdot, :dot, :solid]

    ax_top = Axis(fig[1, 1],
        title="Spatial Resolution - Relative Energy Change",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    for (i, (df, Nx, Nt)) in enumerate(spatial_data)
        t, E_total, _, _ = compute_energy_conservation(df)
        isnothing(t) && continue
        length(E_total) < 2 && continue
        E_total[1] == 0 && continue

        dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
        label = spatial_label_fn(Nx, Nt)
        lines!(ax_top, t, dE_E, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    hlines!(ax_top, [0], color=:black, linestyle=:dash)
    axislegend(ax_top, position=:rt, framevisible=true)

    ax_bottom = Axis(fig[2, 1],
        title="Temporal Resolution - Relative Energy Change",
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false)

    if !isempty(temporal_data)
        Nx_fixed = temporal_data[1][2]
        ax_bottom.title = "Temporal Resolution - Relative Energy Change (Nx=$Nx_fixed)"
    end

    default_temporal_label_fn = let nt_ref=temporal_nt_ref
        if isnothing(nt_ref)
            Nt -> "Nt=$Nt"
        else
            function(nt::Int)
                ratio = nt_ref / nt
                if abs(ratio - 1) < 1e-6
                    "Δt ×1.0"
                else
                    "Δt ×$(round(ratio, digits=1))"
                end
            end
        end
    end
    actual_temporal_label_fn = isnothing(temporal_label_fn) ? default_temporal_label_fn : temporal_label_fn

    for (i, (df, Nx, Nt)) in enumerate(temporal_data)
        t, E_total, _, _ = compute_energy_conservation(df)
        isnothing(t) && continue
        length(E_total) < 2 && continue
        E_total[1] == 0 && continue

        dE_E = (E_total .- E_total[1]) ./ abs(E_total[1])
        label = actual_temporal_label_fn(Nt)
        lines!(ax_bottom, t, dE_E, color=colors[min(i, length(colors))],
               linewidth=2, linestyle=linestyles[min(i, length(linestyles))], label=label)
    end
    hlines!(ax_bottom, [0], color=:black, linestyle=:dash)
    axislegend(ax_bottom, position=:rt, framevisible=true)

    save_path = joinpath(figures_dir, "$save_prefix.png")
    Makie.save(save_path, fig, px_per_unit=2)
    @info "Saved: $save_path"

    return fig
end

let
Nx_files = [
    # "H:/202601_UFDs/Crater_2_2026-04-07_Nx(64)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(64), Xmax(3.85), Nt(209), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-09_Nx(96)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(96), Xmax(3.85), Nt(219), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-07_Nx(128)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(128), Xmax(3.85), Nt(224), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-09_Nx(192)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(192), Xmax(3.85), Nt(229), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-07_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(232), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-09_Nx(384)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(384), Xmax(3.85), Nt(377), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-09_Nx(512)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(512), Xmax(3.85), Nt(671), Tmax(0.09), DM_m(1.00) - Prop.csv",

    "H:/202601_UFDs/Crater_2_2026-04-11_Nx(64)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.0)_rot(0.9), Nx(64), Xmax(3.85), Nt(209), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-11_Nx(96)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.0)_rot(0.9), Nx(96), Xmax(3.85), Nt(219), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-11_Nx(128)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.0)_rot(0.9), Nx(128), Xmax(3.85), Nt(224), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-11_Nx(192)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.0)_rot(0.9), Nx(192), Xmax(3.85), Nt(229), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-11_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.0)_rot(0.9), Nx(256), Xmax(3.85), Nt(232), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/XrayETGs/Crater_2_2026-04-11_Nx(384)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.0)_rot(0.9), Nx(384), Xmax(3.85), Nt(377), Tmax(0.09), DM_m(1.00) - Prop.csv",

    # "H:/202601_UFDs/Crater_2_2026-04-11_Nx(64)_ma(200.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_bs(4)_ac(0)_conv/Crater_2_ma(200.0)_vel(0.0)_rot(0.9), Nx(64), Xmax(7.71), Nt(833), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-11_Nx(96)_ma(200.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_bs(4)_ac(0)_conv/Crater_2_ma(200.0)_vel(0.0)_rot(0.9), Nx(96), Xmax(7.71), Nt(875), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-11_Nx(128)_ma(200.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_bs(4)_ac(0)_conv/Crater_2_ma(200.0)_vel(0.0)_rot(0.9), Nx(128), Xmax(7.71), Nt(896), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/XrayETGs/Crater_2_2026-04-11_Nx(256)_ma(200.0)_mr(1.0)_rr(1.0)_vel(0.0)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_bs(4)_ac(0)_conv/Crater_2_ma(200.0)_vel(0.0)_rot(0.9), Nx(256), Xmax(7.71), Nt(925), Tmax(0.09), DM_m(1.00) - Prop.csv",
]
    
Nt_files = [
    # "H:/202601_UFDs/Crater_2_2026-04-09_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(15), Tmax(0.09), DM_m(1.00) - Prop.csv",
    # "H:/202601_UFDs/Crater_2_2026-04-09_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(29), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-07_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(58), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-07_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(116), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-07_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(232), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-07_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(463), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-07_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(925), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-09_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(1156), Tmax(0.09), DM_m(1.00) - Prop.csv",
    "H:/202601_UFDs/Crater_2_2026-04-09_Nx(256)_ma(50.0)_mr(1.0)_rr(1.0)_vel(0.7)_rot(0.9)_L(40 kpc)_lb(6.0 Gyr)_MW(false)_LMC(false)_bs(4)_ac(0)_conv/Crater_2_ma(50.0)_vel(0.7)_rot(0.9), Nx(256), Xmax(3.85), Nt(2311), Tmax(0.09), DM_m(1.00) - Prop.csv",
]

plot_mass_convergence(Nx_files, Nt_files; temporal_nt_ref = 232)
plot_energy_convergence(Nx_files, Nt_files; temporal_nt_ref = 232)
plot_energy_absolute(Nx_files, Nt_files; temporal_nt_ref = 232)
plot_energy_relative(Nx_files, Nt_files; temporal_nt_ref = 232)
end