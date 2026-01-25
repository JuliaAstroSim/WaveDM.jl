function plot_RC_RAR(dfAcc;
    model = :MW,
    size = (1200, 600),
    Galaxy_id = 0,
    average = false,
    section = 100,
    FDM_mass_ratio = 1.0,
    best_fit_halo_mass = false,
    fontsize = 18,
)
    
    if average
        vc = ustrip.(u"km/s", sqrt.(dfAcc.a_all_averaged * u"m/s^2" .* dfAcc.r * u"kpc"))
    else
        vc = ustrip.(u"km/s", sqrt.(dfAcc.a_all * u"m/s^2" .* dfAcc.r * u"kpc"))
    end
    r_mean, vc_mean, r_std, vc_std = distribution(dfAcc.r, vc; section, uniform_interval=true)

    vc_MOND = ustrip.(u"km/s", sqrt.(dfAcc.a_mond * u"m/s^2" .* dfAcc.r * u"kpc"))
    r_mean, vc_MOND_mean, r_std, vc_MOND_std = distribution(dfAcc.r, vc_MOND; section, uniform_interval=true)

    mass_ratio = mass_fix_ratio * FDM_mass_ratio

    if model == :MW
        dfEilers2019 = load_MW_RC_Eilers2019()
        # dfMroz2019 = load_MW_RC_Mroz2019()
        dfDS_W21 = load_MW_RC_DS_W21()
        dfRCstddev = load_MW_RC_stddev_W21()

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
    
        @info "chi2RC = $(chi2RC)"
        @info "chi2RC_min = $(chi2RC_min)"

        label_scatter = best_fit_halo_mass ? "WaveDM minimum χ²: $(@sprintf("%.2f", chi2RC_min)), mass_ratio: $(@sprintf("%.2f", mass_ratio))" : "WaveDM χ²: $(@sprintf("%.2f", chi2RC))"
        fit_error_min = chi2RC_min
    elseif model == :SPARC_LTGs
        dfSPARC = load_SPARC_LTGs_data()
        df_SPARC_RC = load_SPARC_LTGs_RC()
        df_galaxy_RC = filter(:Galaxy => x->x==dfSPARC.Galaxy[Galaxy_id], df_SPARC_RC)

        chi2RC = chi2reduced(df_galaxy_RC.Vobs, df_galaxy_RC.e_Vobs, df_galaxy_RC.R, vc_mean, r_mean; DOF = 2)

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
    
        @info "chi2RC = $(chi2RC)"
        @info "chi2RC_min = $(chi2RC_min)"

        label_scatter = best_fit_halo_mass ? "WaveDM minimum χ²: $(@sprintf("%.2f", chi2RC_min)), mass_ratio: $(@sprintf("%.2f", mass_ratio))" : "WaveDM χ²: $(@sprintf("%.2f", chi2RC))"
        fit_error_min = chi2RC_min
    else # No observational RC provided, use RAR RC
        MSE_RC = mse(vc_MOND, vc_mean)
        if best_fit_halo_mass
            # Find the minimum MSE by tunning the DM halo mass
            function merit_function(params)
                mass_fix_ratio = params[1]
                merit_chi2RC = mse(vc_MOND, vc_mean * mass_fix_ratio)
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
    
            MSE_RC_min = Optim.minimum(result)[1]
            mass_fix_ratio = Optim.minimizer(result)[1]
        else
            MSE_RC_min = MSE_RC
            mass_fix_ratio = 1
        end

        @info "MSE_RC = $(MSE_RC)"
        @info "MSE_RC_min = $(MSE_RC_min)"

        label_scatter = best_fit_halo_mass ? "WaveDM minimum MSE to RAR RC: $(@sprintf("%.2f", MSE_RC_min)), mass_ratio: $(@sprintf("%.2f", mass_ratio))" : "WaveDM MSE to RAR RC: $(@sprintf("%.2f", MSE_RC))"
        RC_fit_error_min = MSE_RC_min
    end
    

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

    plots = []
    labels = []

    if model == :MW
        l1 = Makie.lines!(ax, dfEilers2019.r, dfEilers2019.v, color = :orange)
        b1 = Makie.band!(ax, dfEilers2019.r, dfEilers2019.v - dfEilers2019.σ_low, dfEilers2019.v + dfEilers2019.σ_high, color = (:orange, 0.2))
        
        s3 = Makie.scatter!(ax, dfRCstddev.r, dfRCstddev.v, marker = :x, color = :gray)
        e3 = Makie.errorbars!(ax, dfRCstddev.r, dfRCstddev.v, dfRCstddev.σ_v, color = :gray)

        l3 = Makie.lines!(dfDS_W21.r, dfDS_W21.v_CDM)

        push!(plots, [l1, b1]); push!(labels, "Observation (Eilers 2019)")
        push!(plots, [s3, e3]); push!(labels, "CDM")
        push!(plots, l3); push!(labels, "Observation (averaged)")

        # s3 = Makie.scatter!(ax, dfMroz2019.r, dfMroz2019.v, color = :blue, markersize = 5.0) # "Mroz 2019"
        # l2 = Makie.lines!(ax, dfRCstddev.r, dfRCstddev.v) # "stddev"
    end

    s1 = Makie.scatter!(ax, r_mean, vc_mean * mass_fix_ratio, color = :red, markersize = 5.0, marker = :x)
    e1 = Makie.errorbars!(ax, r_mean, vc_mean * mass_fix_ratio, vc_std, color = :red)

    s2 = Makie.scatter!(ax, r_mean, vc_MOND_mean, color = :blue, markersize = 5.0)
    e2 = Makie.errorbars!(ax, r_mean, vc_MOND_mean, vc_MOND_std, color = :blue)

    push!(plots, [s1, e1]); push!(labels, )
    push!(plots, [s2, e2]); push!(labels, "RAR")


    Legend(fig[1,1], plots, labels;
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :bottom,
        margin = (10, 10, 10, 10),
    )

    fig, RC_fit_error_min, mass_fix_ratio
end