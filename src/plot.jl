function plot_MW_RC_SPE(dfAcc;
    fileEilers2019 = joinpath(@__DIR__, "../data/MilkyWay/MW_RC_Eilers2019.dat"),
    fileMroz2019 = joinpath(@__DIR__, "../data/MilkyWay/MW_RC_Mroz2019.dat"),
    fileRCstddev = joinpath(@__DIR__, "../data/MilkyWay/MW_RC_stddev_W21.txt"),
    fileDS_W21 = joinpath(@__DIR__, "../data/MilkyWay/MW_dsvel-W21-RC.txt"),
    size = (1200, 600),
    average = false,
    section = 100,
    FDM_mass_ratio = 1.0,
    best_fit_halo_mass = false,
    fontsize = 18,
)
    dfEilers2019 = DataFrame(CSV.File(fileEilers2019; header = false, delim=" ", ignorerepeated = true))
    rename!(dfEilers2019, [:r, :v, :σ_low, :σ_high])

    dfMroz2019 = DataFrame(CSV.File(fileMroz2019; header = false, delim=" ", ignorerepeated = true))
    rename!(dfMroz2019, [:r, :v, :σ_low, :σ_high])

    dfRCstddev = DataFrame(CSV.File(fileRCstddev; header = false, delim=" ", ignorerepeated = true))
    rename!(dfRCstddev, [:id, :r, :v, :σ_v, :σ_r])

    dfDS_W21 = DataFrame(CSV.File(fileDS_W21; header = false, skipto=3, delim="\t", ignorerepeated = true))
    rename!(dfDS_W21, [:r, :v_b, :v_CDM, :v_QUMOND, :v_MOG])

    if average
        vc = ustrip.(u"km/s", sqrt.(dfAcc.a_all_averaged * u"m/s^2" .* dfAcc.r * u"kpc"))
    else
        vc = ustrip.(u"km/s", sqrt.(dfAcc.a_all * u"m/s^2" .* dfAcc.r * u"kpc"))
    end

    r_mean, vc_mean, r_std, vc_std = distribution(dfAcc.r, vc; section, uniform_interval=true)
    chi2RC = chi2reduced(dfRCstddev.v, dfRCstddev.σ_v, dfRCstddev.r, vc_mean, r_mean; DOF = 2)

    if best_fit_halo_mass
        # Find the minimum chi2RC by tunning the DM halo mass

        function merit_function(params)
            mass_fix_ratio = params[1]
            merit_chi2RC = chi2reduced(dfRCstddev.v, dfRCstddev.σ_v, dfRCstddev.r, vc_mean * mass_fix_ratio, r_mean; DOF = 2)
            return merit_chi2RC
        end

        result = Optim.optimize(
            merit_function,
            [0.1], # lower
            [10.0], # upper
            [1.0], # ic
            Optim.Fminbox(),
            Optim.Options(
                # store_trace = true,
                # iterations = 100,
                # outer_iterations = 100,
                # x_tol = 1e-10,
                # outer_x_tol = 1e-10,
            )
        )

        chi2RC_min = Optim.minimum(result)[1]
        mass_fix_ratio = Optim.minimizer(result)[1]
    else
        chi2RC_min = chi2RC
        mass_fix_ratio = 1
    end

    mass_ratio = mass_fix_ratio * FDM_mass_ratio

    @info "chi2RC = $(chi2RC)"
    @info "chi2RC_min = $(chi2RC_min)"


    vc_MOND = ustrip.(u"km/s", sqrt.(dfAcc.a_mond * u"m/s^2" .* dfAcc.r * u"kpc"))
    r_mean, vc_MOND_mean, r_std, vc_MOND_std = distribution(dfAcc.r, vc_MOND; section, uniform_interval=true)

    fig = Figure(; size, fontsize)
    ax = Axis(fig[1,1];
        xlabel = "R [kpc]",
        ylabel = L"v_c [km/s]",
        xminorticksvisible = true,
		# xminorgridvisible = true,
		xminorticks = IntervalsBetween(10),
		yminorticksvisible = true,
		# yminorgridvisible = true,
		yminorticks = IntervalsBetween(10),
    )
    Makie.xlims!(ax, 0, maximum(dfEilers2019.r) * 1.1)
    Makie.ylims!(ax, 0, 300)

    l1 = Makie.lines!(ax, dfEilers2019.r, dfEilers2019.v, color = :orange)
    b1 = Makie.band!(ax, dfEilers2019.r, dfEilers2019.v - dfEilers2019.σ_low, dfEilers2019.v + dfEilers2019.σ_high, color = (:orange, 0.2))
    
    s3 = Makie.scatter!(ax, dfRCstddev.r, dfRCstddev.v, marker = :x, color = :gray)
    e3 = Makie.errorbars!(ax, dfRCstddev.r, dfRCstddev.v, dfRCstddev.σ_v, color = :gray)

    s1 = Makie.scatter!(ax, r_mean, vc_mean * mass_fix_ratio, color = :red, markersize = 5.0, marker = :x)
    e1 = Makie.errorbars!(ax, r_mean, vc_mean * mass_fix_ratio, vc_std, color = :red)

    s2 = Makie.scatter!(ax, r_mean, vc_MOND_mean, color = :blue, markersize = 5.0)
    e2 = Makie.errorbars!(ax, r_mean, vc_MOND_mean, vc_MOND_std, color = :blue)


    # s3 = Makie.scatter!(ax, dfMroz2019.r, dfMroz2019.v, color = :blue, markersize = 5.0)
    # l2 = Makie.lines!(ax, dfRCstddev.r, dfRCstddev.v)

    l3 = Makie.lines!(dfDS_W21.r, dfDS_W21.v_CDM)

    Legend(fig[1,1],
        [
            [s1, e1],
            [s2, e2],
            [s3, e3],
            [l1, b1],
            # s3,
            # l2,
            l3,
        ],
        [
            best_fit_halo_mass ? "WaveDM minimum χ²=$(@sprintf("%.2f", chi2RC_min)) with mass_ratio=$(@sprintf("%.2f", mass_ratio))" : "WaveDM χ²=$(@sprintf("%.2f", chi2RC))",
            "RAR",
            "CDM",
            "Observation (Eilers 2019)",
            # "Mroz 2019",
            # "stddev",
            "Observation (averaged)",
        ];
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :bottom,
        margin = (10, 10, 10, 10),
    )

    fig, chi2RC_min, mass_fix_ratio
end