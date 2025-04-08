function plot_MW_RC_SPE(dfAcc;
    fileEilers2019 = joinpath(@__DIR__, "../data/MilkyWay/MW_RC_Eilers2019.dat"),
    fileMroz2019 = joinpath(@__DIR__, "../data/MilkyWay/MW_RC_Mroz2019.dat"),
    fileRCstddev = joinpath(@__DIR__, "../data/MilkyWay/MW_RC_stddev_W21.txt"),
    fileDS_W21 = joinpath(@__DIR__, "../data/MilkyWay/MW_dsvel-W21-RC.txt"),
    size = (1200, 600),
    average = false,
)
    dfEilers2019 = DataFrame(CSV.File(fileEilers2019; header = false, delim=" ", ignorerepeated = true))
    rename!(dfEilers2019, [:r, :v, :σ_low, :σ_high])

    dfMroz2019 = DataFrame(CSV.File(fileMroz2019; header = false, delim=" ", ignorerepeated = true))
    rename!(dfMroz2019, [:r, :v, :σ_low, :σ_high])

    dfRCstddev = DataFrame(CSV.File(fileRCstddev; header = false, delim=" ", ignorerepeated = true))
    rename!(dfRCstddev, [:id, :r, :v, :σ_v, :σ_r])

    dfDS_W21 = DataFrame(CSV.File(fileDS_W21; header = false, skipto=3, delim="\t", ignorerepeated = true))
    rename!(dfDS_W21, [:r, :v_b, :v_CDM, :v_QUMOND, :v_MOG])

    fig = Figure(; size)
    ax = Axis(fig[1,1];
        xlabel = "r [kpc]",
        ylabel = "v_rot [km/s]",
    )
    Makie.xlims!(ax, 0, maximum(dfEilers2019.r) * 1.1)
    Makie.ylims!(ax, 0, 300)

    if average
        vc = ustrip.(u"km/s", sqrt.(dfAcc.a_all_averaged * u"m/s^2" .* dfAcc.r * u"kpc"))
    else
        vc = ustrip.(u"km/s", sqrt.(dfAcc.a_all * u"m/s^2" .* dfAcc.r * u"kpc"))
    end
    s1 = Makie.scatter!(ax, dfAcc.r, vc, color = :red, markersize = 5.0, marker = :x)

    vc_MOND = ustrip.(u"km/s", sqrt.(dfAcc.a_mond * u"m/s^2" .* dfAcc.r * u"kpc"))
    s2 = Makie.scatter!(ax, dfAcc.r, vc_MOND, color = :blue, markersize = 5.0)
    
    l1 = Makie.lines!(ax, dfEilers2019.r, dfEilers2019.v, color = :black)
    b1 = Makie.band!(ax, dfEilers2019.r, dfEilers2019.v - dfEilers2019.σ_low, dfEilers2019.v + dfEilers2019.σ_high, color = (:black, 0.2))

    # reduced chi2
    # chi2RC = chi2reduced(dfEilers2019.v, dfEilers2019.σ_high - dfEilers2019.σ_low, dfEilers2019.r, vc, dfAcc.r; DOF = 2)
    # Use the averaged dispersion data
    chi2RC = chi2reduced(dfRCstddev.v, dfRCstddev.σ_v, dfRCstddev.r, vc, dfAcc.r; DOF = 2)

    s3 = Makie.scatter!(ax, dfRCstddev.r, dfRCstddev.v, marker = :x, color = :black)
    e3 = Makie.errorbars!(ax, dfRCstddev.r, dfRCstddev.v, dfRCstddev.σ_v, color = :black)

    # s3 = Makie.scatter!(ax, dfMroz2019.r, dfMroz2019.v, color = :blue, markersize = 5.0)
    # l2 = Makie.lines!(ax, dfRCstddev.r, dfRCstddev.v)

    l3 = Makie.lines!(dfDS_W21.r, dfDS_W21.v_CDM)

    Legend(fig[1,1],
        [
            s1,
            s2,
            [l1, b1],
            # s3,
            # l2,
            l3,
            [s3, e3],
        ],
        [
            "WaveDM χ²=$(@sprintf("%.2f", chi2RC))",
            "RAR",
            "Eilers 2019",
            # "Mroz 2019",
            # "stddev",
            "CDM",
            "Averaged",
        ];
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :bottom,
        margin = (10, 10, 10, 10),
    )

    fig, chi2RC
end