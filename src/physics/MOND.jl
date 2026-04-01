RAR(a_b, a0) = a_b / (1 - exp(-sqrt(abs(a_b)/a0)))

Milgrom(::Missing, ::Any, ::Any) = missing
Milgrom(::Any, ::Missing, ::Any) = missing
Milgrom(a_b, a0, nuIndex) = nu(abs(a_b)/a0, nuIndex) * a_b

RAR_inv(a_b, a_MOND) = a_MOND > a_b ? a_b / (log(1-a_b/a_MOND))^2 : missing
Milgrom_inv(a_b, a_MOND, nuIndex) = a_MOND > a_b ? a_b * (((2*(a_MOND/a_b)^nuIndex - 1)^2 - 1) / 4)^(1/nuIndex) : missing

# Lsq fitting models
# modelRAR(x, p) = x ./ (1 .- exp.(-sqrt.(abs.(x)/abs(p[1]))))
modelRAR(x, p) = x ./ (1 .- exp.(-sqrt.(abs.(x)/p[1])))


function compute_RAR(dfAcc;
    minR = 8.0u"kpc",
    maxR = 15.0u"kpc",
    plotMinR = 0.0u"kpc",
    plotMaxR = 30.0u"kpc",
    nuIndex = 2.0,
    zoom_a_max = 2e-10,
    average = false,
    mass_fix_ratio = 1.0,
    a0only = false,
    RARonly = false,
    size = (1920, 1080),
)
    @show minR
    @show maxR
    @show plotMinR
    @show plotMaxR

    indices = plotMinR .< (dfAcc.r * u"kpc") .< plotMaxR
    r = dfAcc.r[indices]
    a_b = dfAcc.a_b[indices]
    if average
        a_all = dfAcc.a_all_averaged[indices] * mass_fix_ratio
    else
        a_all = dfAcc.a_all[indices] * mass_fix_ratio
    end
    a_mond = dfAcc.a_mond[indices]

    # interested region where there might have obseravtional data
    indices_data = minR .< (dfAcc.r * u"kpc") .< maxR
    r_data = dfAcc.r[indices_data]
    a_b_data = dfAcc.a_b[indices_data]
    if average
        a_all_data = dfAcc.a_all_averaged[indices_data] * mass_fix_ratio
    else
        a_all_data = dfAcc.a_all[indices_data] * mass_fix_ratio
    end
    a_mond_data = dfAcc.a_mond[indices_data]


    a0_RAR = RAR_inv.(a_b, a_all)
    a0_Milgrom = Milgrom_inv.(a_b, a_all, nuIndex)

    a0_RAR_data = RAR_inv.(a_b_data, a_all_data)
    a0_Milgrom_data = Milgrom_inv.(a_b_data, a_all_data, nuIndex)

    a0_RAR_data_mean = mean(skipmissing(a0_RAR_data))
    a0_Milgrom_data_mean = mean(skipmissing(a0_Milgrom_data))

    a0_H0 = ustrip(u"m/s^2", C.c * C.H / (2π)) # 1.1447167528358282e-10 m s^-2

    a_mond_RAR = RAR.(a_b, a0_RAR_data_mean)
    a_mond_Milgrom = Milgrom.(a_b, a0_Milgrom_data_mean, nuIndex)
    #TODO: RAR 用 Milgrom a_0 算，另一个同理？
    a_mond_RAR_H0 = RAR.(a_b, a0_H0)
    a_mond_Milgrom_H0 = Milgrom.(a_b, a0_H0, nuIndex)



    f = Figure(; size);

    axACC0 = Axis(f[1,1];
        xlabel = "r [kpc]",
        ylabel = "a₀ [m/s²]",
        yscale = log10,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
    )
    Makie.xlims!(axACC0, ustrip(u"kpc", plotMinR), ustrip(u"kpc", plotMaxR))
    s1 = Makie.scatter!(axACC0, r, a0_RAR, color=:red, markersize = 5.0)
    hl1 = Makie.hlines!(axACC0, a0_RAR_data_mean, color=:red)

    if !RARonly
        # MOND
        s2 = Makie.scatter!(axACC0, r, a0_Milgrom, color=:orange, markersize = 5.0)
        hl2 = Makie.hlines!(axACC0, a0_Milgrom_data_mean, color=:orange)
    end
    # hl3 = Makie.hlines!(axACC0, a0_H0, color=:brown)
    hl4 = Makie.hlines!(axACC0, ustrip(u"m/s^2", C.ACC0), color=:blue)

    # Makie.vlines!(axACC0, ustrip(u"kpc", minR), color=(:black, 0.5))
    # Makie.vlines!(axACC0, ustrip(u"kpc", maxR), color=(:black, 0.5))

    Makie.vspan!(axACC0, ustrip(u"kpc", plotMinR), ustrip(u"kpc", minR), color=(:black, 0.1))
    Makie.vspan!(axACC0, ustrip(u"kpc", maxR), ustrip(u"kpc", plotMaxR), color=(:black, 0.1))

    if RARonly
        Legend(f[1,1],
            [s1, hl1, hl4],
            [
                "RAR",
                # "Milgrom 1983",
                "mean(RAR) = $(@sprintf("%.4e m/s²", a0_RAR_data_mean))",
                # "mean(Milgrom 1983) = $(@sprintf("%.4e m/s²", a0_Milgrom_data_mean))",
                # "Hubble a₀ ≈ cH₀/(2π) = $(a0_H0)",
                "a₀ = 1.2e-10 m/s²",
            ],
            tellheight = false,
            tellwidth = false,
            halign = :right,
            valign = :top,
            margin = (10, 10, 10, 10),
        )
    else
        #MOND
        Legend(f[1,1],
            [s1, s2, hl1, hl2, hl4],
            [
                "RAR",
                "Milgrom 1983",
                "mean(RAR) = $(@sprintf("%.4e m/s²", a0_RAR_data_mean))",
                "mean(Milgrom 1983) = $(@sprintf("%.4e m/s²", a0_Milgrom_data_mean))",
                # "Hubble a₀ ≈ cH₀/(2π) = $(a0_H0)",
                "a₀ = 1.2e-10 m/s²",
            ],
            tellheight = false,
            tellwidth = false,
            halign = :right,
            valign = :top,
            margin = (10, 10, 10, 10),
        )
    end

    if a0only
        return f
    end

    axACC0_zoom = Axis(f[1,2];
        xlabel = "r [kpc]",
        ylabel = "a₀ [m/s²]",
        xscale = log10,
        yscale = log10,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
    )
    Makie.xlims!(axACC0_zoom, ustrip(u"kpc", minR), ustrip(u"kpc", maxR))
    Makie.ylims!(axACC0_zoom, min(minimum(skipmissing(a0_RAR_data)), minimum(skipmissing(a0_Milgrom_data))), max(maximum(skipmissing(a0_RAR_data)), maximum(skipmissing(a0_Milgrom_data))))
    s1 = Makie.scatter!(axACC0_zoom, r, a0_RAR, color=:blue, markersize = 5.0)
    s2 = Makie.scatter!(axACC0_zoom, r, a0_Milgrom, color=:red, markersize = 5.0)
    hl1 = Makie.hlines!(axACC0_zoom, a0_RAR_data_mean, color=:blue)
    hl2 = Makie.hlines!(axACC0_zoom, a0_Milgrom_data_mean, color=:red)
    hl3 = Makie.hlines!(axACC0_zoom, a0_H0, color=:brown)
    hl4 = Makie.hlines!(axACC0_zoom, ustrip(u"m/s^2", C.ACC0), color=:orange)



    axModels = Axis(f[2,1];
        xlabel = "r [kpc]",
        ylabel = "acc [m/s²]",
        yscale = log10,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
    )
    Makie.xlims!(axModels, ustrip(u"kpc", plotMinR), ustrip(u"kpc", plotMaxR))
    s1 = Makie.scatter!(axModels, r, a_b, color = :black, markersize = 5.0)
    s2 = Makie.scatter!(axModels, r, a_all, color = :green, markersize = 5.0)
    s3 = Makie.scatter!(axModels, r, a_mond, color = :orange, markersize = 5.0)
    # s4 = Makie.scatter!(axModels, r, a_mond_RAR, color = :blue, markersize = 5.0)
    s5 = Makie.scatter!(axModels, r, a_mond_Milgrom, color = :red, markersize = 5.0)
    # s6 = Makie.scatter!(axModels, r, a_mond_RAR_H0, color = :brown, markersize = 5.0)
    # s7 = Makie.scatter!(axModels, r, a_mond_Milgrom_H0, color = :cyan, markersize = 5.0)

    hl1 = Makie.hlines!(axModels, ustrip(u"m/s^2", C.ACC0), color=:orange)
    hl2 = Makie.hlines!(axModels, a0_RAR_data_mean, color=:blue)
    hl3 = Makie.hlines!(axModels, a0_Milgrom_data_mean, color=:red)
    hl4 = Makie.hlines!(axModels, a0_H0, color=:brown)

    Makie.vspan!(axModels, ustrip(u"kpc", plotMinR), ustrip(u"kpc", minR), color=(:black, 0.1))
    Makie.vspan!(axModels, ustrip(u"kpc", maxR), ustrip(u"kpc", plotMaxR), color=(:black, 0.1))

    Legend(f[2,1],
        [
            s1,
            s2,
            [s3, hl1],
            # [s4, hl2],
            [s5, hl3],
            # [s6, hl4],
            # s7
        ],
        [
            "baryon",
            "all (baryon + WaveDM)",
            "MOND (a₀ = 1.2e-10 m/s²)",
            # "RAR (a₀ = $(@sprintf("%.4e m/s²", a0_RAR_data_mean)), RAR)",
            "MOND (a₀ = $(@sprintf("%.4e m/s²", a0_Milgrom_data_mean)), Milgrom 1983)",
            # "RAR (a₀ = $(@sprintf("%.4e m/s²", a0_H0)), H0)",
            # "MOND (a₀ = $(@sprintf("%.4e m/s²", a0_H0)), H0)",
        ],
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
        margin = (20, 20, 20, 20),
        # markersize = 50,
        # labelsize = 30,
    )

    axModels_zoom = Axis(f[2,2];
        xlabel = "r [kpc]",
        ylabel = "acc [m/s²]",
        xscale = log10,
        yscale = log10,
        xminorticksvisible = true,
        xminorgridvisible = true,
        yminorticksvisible = true,
        yminorgridvisible = true,
    )
    Makie.xlims!(axModels_zoom, ustrip(u"kpc", minR), ustrip(u"kpc", maxR))
    Makie.ylims!(axModels_zoom, minimum(a_b_data), zoom_a_max)
    s1 = Makie.scatter!(axModels_zoom, r, a_b, color = :black, markersize = 5.0)
    s2 = Makie.scatter!(axModels_zoom, r, a_all, color = :green, markersize = 5.0)
    s3 = Makie.scatter!(axModels_zoom, r, a_mond, color = :orange, markersize = 5.0)
    # s4 = Makie.scatter!(axModels_zoom, r, a_mond_RAR, color = :blue, markersize = 5.0)
    s5 = Makie.scatter!(axModels_zoom, r, a_mond_Milgrom, color = :red, markersize = 5.0)
    # s6 = Makie.scatter!(axModels_zoom, r, a_mond_RAR_H0, color = :brown, markersize = 5.0)
    # s7 = Makie.scatter!(axModels_zoom, r, a_mond_Milgrom_H0, color = :cyan, markersize = 5.0)

    f
end

function plot_acc_RAR!(ax, dfAcc;
    rMin = 0.0, # kpc
    rMax = 20.0, # kpc
    xLim = (1e-13, 1e-9),
    yLim = (1e-13, 1e-9),
    average = false,
    mass_fix_ratio = 1.0,
    marker = :circle,
    markersize = 5.0,
    color = :red,
    flag_lines = true,
    kw...
)
    indices = rMin .<= dfAcc.r .<= rMax
    if average
        s1 = Makie.scatter!(ax, dfAcc.a_b[indices], dfAcc.a_all_averaged[indices] * mass_fix_ratio; color, markersize, marker, kw...)
    else
        s1 = Makie.scatter!(ax, dfAcc.a_b[indices], dfAcc.a_all[indices] * mass_fix_ratio; color, markersize, marker, kw...)
    end
    if flag_lines
        l1 = Makie.lines!(ax, collect(xLim), collect(xLim); color = :black, linestyle = :dot) # line of unity
    else
        l1 = nothing
    end
    Makie.xlims!(ax, xLim)
    Makie.ylims!(ax, yLim)

    # a0 = 1.2e-10
    x = collect(LinRange(xLim..., 100))
    if flag_lines
        y = RAR.(x, 1.2e-10)
        l2 = Makie.lines!(ax, x, y; color = :blue)
    else
        l2 = nothing
    end

    # Lsq fitting of a0
    if average
        @time "Lsq fitting" fitRAR = curve_fit(modelRAR, dfAcc.a_b[indices] * 1e10, dfAcc.a_all_averaged[indices] * mass_fix_ratio * 1e10, [1.2];
            lower = [1e-2],
            upper = [1e2],
            # show_trace=true,
        )
    else
        @time "Lsq fitting" fitRAR = curve_fit(modelRAR, dfAcc.a_b[indices] * 1e10, dfAcc.a_all[indices] * mass_fix_ratio * 1e10, [1.2];
            lower = [1e-2],
            upper = [1e2],
            # show_trace=true,
        )
    end

    y = RAR.(x, fitRAR.param[1] * 1e-10)
    l3 = Makie.lines!(ax, x, y; color, linestyle = :dash)

    return s1, l1, l2, l3, fitRAR
end

function plot_acc_RAR(dfAcc;
    size = (800, 600),
    rMin = 0.0, # kpc
    rMax = 20.0, # kpc
    xLim = (1e-13, 1e-9),
    yLim = (1e-13, 1e-9),
    title = "",
    markersize = 5.0,
    average = false,
    mass_fix_ratio = 1.0,
    kw...
)
    fig = Figure(; size);
    ax = Axis(fig[1,1];
        xlabel = L"\log_{10}(g_{bar} [m/s^2])",
        ylabel = L"\log_{10}(g_{obs} [m/s^2])",
        xscale = log10,
        yscale = log10,
        title,
    )

    s1, l1, l2, l3, fitRAR = plot_acc_RAR!(ax, dfAcc; rMin, rMax, xLim, yLim, markersize, average, mass_fix_ratio, kw...)

    Legend(fig[1,1],
        [
            s1,
            l1,
            l2,
            l3,
        ],
        [
            "simulation",
            "unity",
            "RAR (a₀ = 1.20e-10 [m/s²])",
            "RAR (a₀ = $(@sprintf("%.2f", fitRAR.param[1]))e-10 [m/s²])",
        ];
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :bottom,
        margin = (10, 10, 10, 10),
    )

    fig
end