"""
$(TYPEDSIGNATURES)

Compute density profile fitting error.
This function encapsulates the profile fitting error computation from SPE3D_MOND.
"""
function compute_profile_fit_error(r_mass_center, rho, length_astro, Δ,
    target_profile_model, target_profile_ρ0, target_profile_rs, target_profile_α,
    target_profile_β, target_profile_γ, density_astro, uniform_interval)
    
    r_filter_3D = 0 .< r_mass_center .<= target_profile_rs * 1.0 / length_astro  # 3D filter
    r_mean, rho_mean, r_std, rho_std = distribution(r_mass_center[r_filter_3D], collect(rho)[r_filter_3D];
        section = floor(Int, target_profile_rs * 1.0 / length_astro / Δ[1]),
        uniform_interval,
    );

    if target_profile_model == :dwarf_gNFW
        model_halo = gNFW(target_profile_β, target_profile_ρ0, target_profile_rs)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    elseif target_profile_model == :dwarf_NFW
        model_halo = NFW(target_profile_ρ0, target_profile_rs)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    elseif target_profile_model == :dwarf_Zhao
        model_halo = Zhao(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    elseif target_profile_model == :dwarf_ZhaoQ #TODO consider Q
        model_halo = Zhao(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ)
        # model_halo = ZhaoQ(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ, target_profile_Q)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    end

    r_filter = 0 .< r_mean .<= target_profile_rs / length_astro  # 1D filter
    current_fit_error = sum((log.(rho_mean[r_filter]) .- log.(ρ_halo[r_filter])).^2) / sum(r_filter) # mse
    
    return current_fit_error
end

"""
$(TYPEDSIGNATURES)

Compute rotation curve fitting error.
This function encapsulates the RC fitting error computation from SPE3D_MOND.
"""
function compute_rc_fit_error(r_mass_center, ax_all, ay_all, az_all, xxx, yyy, zzz,
    rho_max_id, target_profile_rs, length_astro, Δ, velocity_astro, uL,
    df_CO_RC, uniform_interval)
    
    r_filter_3D = 0 .< r_mass_center .<= target_profile_rs * 1.0 / length_astro  # 3D filter
    ar_all = vec_cartesian_to_spherical(ax_all, ay_all, az_all, xxx.-xxx[rho_max_id], yyy.-yyy[rho_max_id], zzz.-zzz[rho_max_id], r_mass_center)[1];
    indices = findall(x->isnan(x), ar_all)
    ar_all[indices] .= 0
    r_mean, ar_all_mean, _, ar_all_std = distribution(r_mass_center[r_filter_3D], abs.(ar_all[r_filter_3D]); section = ceil(Int, target_profile_rs / length_astro / Δ[1]), uniform_interval);
    vel_rot_all_mean = ustrip.(u"km/s", sqrt.(ar_all_mean .* r_mean) * velocity_astro)
    vel_rot_all_u = ustrip.(u"km/s", sqrt.((ar_all_mean + ar_all_std) .* r_mean) * velocity_astro)
    ar_all_d = ar_all_mean - ar_all_std
    ar_all_d[ar_all_d .< 0] .= 0.0
    vel_rot_all_d = ustrip.(u"km/s", sqrt.(ar_all_d .* r_mean) * velocity_astro)
    current_fit_error = chi2reduced(df_CO_RC.vel, df_CO_RC.vel_e, df_CO_RC.r, vel_rot_all_mean, vel_rot_all_u-vel_rot_all_d, r_mean * uL)
    
    return current_fit_error
end

"""
$(TYPEDSIGNATURES)

Compute beta_star (inner slope) from density profile.
This function encapsulates the beta_star computation from SPE3D_MOND.
"""
function compute_beta_star(r_mass_center, rho, target_beta_star_r_min, target_beta_star_r_max, length_astro)
    r_filter = target_beta_star_r_min/length_astro .<= r_mass_center .<= target_beta_star_r_max/length_astro  # 3D filter
    model(t, p) = p[1] .* t .+ p[2]
    log_r = log.(r_mass_center[r_filter])
    log_ρ = log.(collect(rho)[r_filter])
    initial_guess = [-1.0, 0.0]
    fit = curve_fit(model, log_r, log_ρ, initial_guess)
    current_beta_star, intercept_inner = fit.param
    return current_beta_star
end

"""
$(TYPEDSIGNATURES)

Update best fit snapshot based on fitting errors.
This function encapsulates the optimization update logic from SPE3D_MOND.
"""
function update_best_fit!(best_fit_error, current_fit_error, t, i, time_astro,
    best_fit_ψ, ψ, best_fit_ψ_last_t, ψ_last_t, best_fit_Φ_all, Φ_all,
    target_beta_star, beta_star_error_threshold, fig, outputdir, title, suffix,
    r_mass_center, rho, length_astro, target_beta_star_r_min, target_beta_star_r_max)
    
    if isnan(target_beta_star)  # constrain the profile only
        if current_fit_error < best_fit_error
            best_fit_error = current_fit_error
            best_fit_t = t[i] * time_astro
            best_fit_ψ .= ψ
            best_fit_ψ_last_t .= ψ_last_t
            best_fit_Φ_all .= Φ_all

            Makie.save(joinpath(outputdir, "$(title), $(suffix) - Overview Prop best fit.png"), fig)
        end
    else  # constrain both profile and slope
        current_beta_star = compute_beta_star(r_mass_center, rho, target_beta_star_r_min, target_beta_star_r_max, length_astro)
        current_beta_star_fit_error = current_beta_star - target_beta_star
        
        if current_fit_error < best_fit_error && abs(current_beta_star_fit_error) < beta_star_error_threshold
            best_fit_error = current_fit_error
            best_fit_t = t[i] * time_astro
            best_fit_ψ .= ψ
            best_fit_ψ_last_t .= ψ_last_t
            best_fit_Φ_all .= Φ_all

            best_fit_beta_star_error = current_beta_star_fit_error
            best_fit_beta_star = current_beta_star

            Makie.save(joinpath(outputdir, "$(title), $(suffix) - Overview Prop best fit.png"), fig)
        end
    end
    
    return best_fit_error, best_fit_t, best_fit_beta_star_error, best_fit_beta_star
end

# Export functions
export compute_profile_fit_error, compute_rc_fit_error
export update_best_fit!
