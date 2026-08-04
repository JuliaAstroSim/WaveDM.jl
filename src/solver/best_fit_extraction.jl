"""
$(TYPEDSIGNATURES)

Compute density profile fitting error.
This function encapsulates the profile fitting error computation from SPE3D_waveDM.
"""
function compute_profile_fit_error(r_mass_center, rho, length_astro, Δ, density_astro, profile_config::ProfileFitConfig, target_profile_model, uniform_interval)
    target_profile_ρ0 = profile_config.target_profile_ρ0
    target_profile_rs = profile_config.target_profile_rs
    target_profile_α = profile_config.target_profile_α
    target_profile_β = profile_config.target_profile_β
    target_profile_γ = profile_config.target_profile_γ
    
    r_filter_3D = 0 .< r_mass_center .<= target_profile_rs * 1.0 / length_astro  # 3D filter
    r_mean, rho_mean, r_std, rho_std = distribution(r_mass_center[r_filter_3D], collect(rho)[r_filter_3D];
        section = floor(Int, target_profile_rs * 1.0 / length_astro / Δ[1]),
        uniform_interval,
    );

    if target_profile_model == :dwarf_gNFW
        model_halo = gNFW(target_profile_β, target_profile_ρ0, target_profile_rs)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    elseif target_profile_model == :dwarf_NFW || target_profile_model == :NFW
        model_halo = NFW(target_profile_ρ0, target_profile_rs)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    elseif target_profile_model == :dwarf_Zhao || target_profile_model == :Zhao || target_profile_model == :MW
        model_halo = Zhao(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    elseif target_profile_model == :dwarf_ZhaoQ #TODO consider Q
        model_halo = Zhao(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ)
        # model_halo = ZhaoQ(target_profile_ρ0, target_profile_rs, target_profile_α, target_profile_β, target_profile_γ, target_profile_Q)
        ρ_halo = upreferred.(GalacticDynamics.density.(model_halo, r_mean * length_astro) / density_astro)
    else
        error("Unsupported target_profile_model :$(target_profile_model)")
    end

    r_filter = 0 .< r_mean .<= target_profile_rs / length_astro  # 1D filter
    current_fit_error = mse(log.(rho_mean[r_filter]), log.(ρ_halo[r_filter]))
    
    return current_fit_error
end

"""
$(TYPEDSIGNATURES)

Compute rotation curve fitting error.
This function encapsulates the RC fitting error computation from SPE3D_waveDM.
"""
function compute_rc_fit_error(r_mass_center, ax_all, ay_all, az_all, xxx, yyy, zzz, rho_max_id, length_astro, Δ, config_units::AstroUnitsConfig, df_CO_RC, uniform_interval)
    target_profile_rs = 0.1  # TODO: Get this from profile config
    velocity_astro = config_units.velocity_astro
    uL = config_units.uL
    
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
This function encapsulates the beta_star computation from SPE3D_waveDM.
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
This function encapsulates the optimization update logic from SPE3D_waveDM.
"""
function update_best_fit!(best_fit_error, best_fit_t, best_fit_beta_star_error, best_fit_beta_star, current_beta_star, current_fit_error, t, i, time_astro, best_fit_ψ, ψ, best_fit_ψ_last_t, ψ_last_t, best_fit_Φ_all, Φ_all, best_fit_a_all, a_all,
                          rc_config::RCFitConfig, fig, outputdir, title, suffix, r_mass_center, rho, length_astro)
    target_beta_star = rc_config.target_beta_star
    beta_star_error_threshold = rc_config.beta_star_error_threshold
    target_beta_star_r_min = rc_config.target_beta_star_r_min
    target_beta_star_r_max = rc_config.target_beta_star_r_max
    
    if isnan(target_beta_star)  # constrain the profile only
        if current_fit_error < best_fit_error
            best_fit_error = current_fit_error
            best_fit_t = t[i] * time_astro
            best_fit_ψ .= ψ
            best_fit_ψ_last_t .= ψ_last_t
            best_fit_Φ_all .= Φ_all
            best_fit_a_all .= a_all

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
            best_fit_a_all .= a_all

            best_fit_beta_star_error = current_beta_star_fit_error
            best_fit_beta_star = current_beta_star

            Makie.save(joinpath(outputdir, "$(title), $(suffix) - Overview Prop best fit.png"), fig)
        end
    end
    
    return best_fit_error, best_fit_t, best_fit_beta_star_error, best_fit_beta_star, current_beta_star
end


"""
$(TYPEDSIGNATURES)

Compute spherically-averaged radial, zenithal, azimuthal velocity
moments in shells defined by `r_bins` (in code units of `length_astro`).

Returns a NamedTuple with the per-shell mean and standard deviation
of v_r, v_θ, v_φ.

Implementation notes
--------------------
* Velocity is reconstructed from the wavefunction phase via
  `v = ∇arg(ψ)`, *not* from mass-weighted `p = ρ v = Im(ψ* ∇ψ)`,
  because the phase-gradient velocity is the natural Madelung
  variable for a non-interacting condensate (irrotational flow).
* The center is taken as the density maximum `rho_max_id`, which
  matches the convention used elsewhere in the loop.
"""
function compute_velocity_field_shells(
    ψ, xxx, yyy, zzz, r_mass_center, rho_max_id,
    r_bins; # Vector{Float64} of bin edges in code units
)
    # Madelung velocity field: v = ∇S / m = ∇arg(ψ) / m
    phase = angle.(ψ .+ eps(ComplexF64))  # eps() protects against ψ=0

    Nx, Ny, Nz = size(ψ)
    Δx, Δy, Δz = xxx[2,1,1] - xxx[1,1,1], yyy[1,2,1] - yyy[1,1,1], zzz[1,1,2] - zzz[1,1,1]

    # Central differences on the grid (FFT-based would be cleaner; see TODO)
    vx = similar(real.(ψ)); vy = similar(vx); vz = similar(vx)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 2:Nx-1
        vx[i,j,k] = (phase[i+1,j,k] - phase[i-1,j,k]) / (2Δx)
    end
    @inbounds for k in 1:Nz, j in 2:Ny-1, i in 1:Nx
        vy[i,j,k] = (phase[i,j+1,k] - phase[i,j-1,k]) / (2Δy)
    end
    @inbounds for k in 2:Nz-1, j in 1:Ny, i in 1:Nx
        vz[i,j,k] = (phase[i,j,k+1] - phase[i,j,k-1]) / (2Δz)
    end

    # Spherical velocity components about rho_max center
    x_rel = xxx .- xxx[rho_max_id]
    y_rel = yyy .- yyy[rho_max_id]
    z_rel = zzz .- zzz[rho_max_id]
    r_sph = sqrt.(x_rel.^2 .+ y_rel.^2 .+ z_rel.^2)
    eps_r = eps(Float64)
    r_sph_safe = max.(r_sph, eps_r)

    # Unit vectors
    sr_x = x_rel ./ r_sph_safe
    sr_y = y_rel ./ r_sph_safe
    sr_z = z_rel ./ r_sph_safe
    # θ̂ = (cosθ cosφ, cosθ sinφ, -sinθ) in code coords (z is polar axis):
    # θ̂ = (z·x/r², z·y/r², -ρ²_xy/r²)
    theta_x = z_rel .* x_rel ./ r_sph_safe.^2
    theta_y = z_rel .* y_rel ./ r_sph_safe.^2
    theta_z = -(x_rel.^2 .+ y_rel.^2) ./ r_sph_safe.^2

    # v_r = v · r̂
    v_r = sr_x .* vx .+ sr_y .* vy .+ sr_z .* vz
    # v_θ = v · θ̂
    v_θ = theta_x .* vx .+ theta_y .* vy .+ theta_z .* vz
    # v_φ = v · φ̂, with φ̂ = (-y, x, 0)/ρ_xy
    rho_xy = sqrt.(x_rel.^2 .+ y_rel.^2)
    rho_xy_safe = max.(rho_xy, eps_r)
    phi_x = -y_rel ./ rho_xy_safe
    phi_y = x_rel ./ rho_xy_safe
    phi_z = zero(eltype(rho_xy_safe))
    v_φ = phi_x .* vx .+ phi_y .* vy .+ phi_z .* vz

    # Bin by r (r_bins in code units; r_sph also in code units)
    section_count = length(r_bins) - 1
    r_flat  = vec(r_sph)
    vr_flat = vec(v_r)
    vθ_flat = vec(v_θ)
    vφ_flat = vec(v_φ)

    r_mean, vr_mean, _, vr_std = distribution(r_flat, vr_flat;
        section = section_count, uniform_interval = true)
    _, vθ_mean, _, vθ_std = distribution(r_flat, vθ_flat;
        section = section_count, uniform_interval = true)
    _, vφ_mean, _, vφ_std = distribution(r_flat, vφ_flat;
        section = section_count, uniform_interval = true)

    return (
        r_mean = r_mean,
        vr_mean = vr_mean, vr_std = vr_std,
        vθ_mean = vθ_mean, vθ_std = vθ_std,
        vφ_mean = vφ_mean, vφ_std = vφ_std,
    )
end

"""
$(TYPEDSIGNATURES)

Compute the mass-weighted center of the halo, optionally restricted to
bound cells (ρ > `rho_bound_threshold`). This is the V4 "bound mass center"
diagnostic.

Returns `(xc, yc, zc)` in code units.
"""
function compute_mass_weighted_center(
    rho, xxx, yyy, zzz; rho_bound_threshold = 0.0)
    mask = rho .> rho_bound_threshold
    sum_rho = sum(rho[mask])
    if sum_rho == 0
        _, id = findmax(rho)
        return (xxx[id], yyy[id], zzz[id])
    end
    xc = sum(xxx[mask] .* rho[mask]) / sum_rho
    yc = sum(yyy[mask] .* rho[mask]) / sum_rho
    zc = sum(zzz[mask] .* rho[mask]) / sum_rho
    return (xc, yc, zc)
end

"""
$(TYPEDSIGNATURES)

Compute interior (within `r_interior_kpc` kpc) mass. This is the V2
"interior mass vs. time" diagnostic.

Returns interior mass in solar masses.
"""
function compute_interior_mass(rho, r_mass_center, unit_cell_volume, mass_astro;
                                r_interior_kpc = 1.0, uL = 1.0)
    r_max_code = r_interior_kpc / uL
    mask = r_mass_center .<= r_max_code
    return sum(rho[mask]) * unit_cell_volume * mass_astro
end

"""
$(TYPEDSIGNATURES)

Per-timestep diagnostic bundle, returned by `compute_diagnostic_snapshot`.

The fields cover every quantity needed for the cuspy-window / β*(t) /
velocity-field-modulation analysis (Modules A, B, C, D in §10.2 of
the paper planning note).
"""
struct DiagnosticSnapshot
    t::Float64                          # code-unit time
    rho_c::Float64                      # peak density
    xc_mass::Float64                    # mass-weighted center x (code units)
    yc_mass::Float64
    zc_mass::Float64
    bound_mass::Float64                 # bound halo mass (solar masses)
    interior_mass::Float64              # interior mass within r_interior_kpc
    beta_star_rhomax::Float64           # β* about ρ_max
    beta_star_masscenter::Float64       # β* about mass-weighted center (V4)
    beta_star_rhomax_inner::Float64     # β* in tighter interval [0.3, 0.5] kpc
    virial_V::Float64
    virial_Q::Float64
    virial_K::Float64
end

"""
$(TYPEDSIGNATURES)

Compute one diagnostic snapshot from the current SPE state.

This is called once per timestep (after `extract_min_t`) to accumulate
the time series β*(t), ρ_c(t), mass-center coordinates, etc.
"""
function compute_diagnostic_snapshot(
    ψ, rho, xxx, yyy, zzz, r_mass_center, rho_max_id,
    t_current_code,                    # current SPE time in code units
    unit_cell_volume, mass_astro, length_astro,
    QE, KE, PE_abs,
    r_bins;                              # velocity shell bin edges (code units)
    rho_bound_threshold = 0.0,
    r_interior_kpc = 1.0,
    beta_star_r_min_kpc = 0.3,
    beta_star_r_max_kpc = 0.8,
)
    rho_c = rho[rho_max_id]
    xc_m, yc_m, zc_m = compute_mass_weighted_center(rho, xxx, yyy, zzz;
        rho_bound_threshold = rho_bound_threshold)

    rho_thresh = max(rho_bound_threshold, rho_c * 1e-3)
    bound_mass = sum(rho[rho .> rho_thresh]) * unit_cell_volume * mass_astro

    interior_mass = compute_interior_mass(rho, r_mass_center,
        unit_cell_volume, mass_astro;
        r_interior_kpc = r_interior_kpc, uL = length_astro)

    beta_star_rm = compute_beta_star(r_mass_center, rho,
        beta_star_r_min_kpc, beta_star_r_max_kpc, length_astro)

    # β* about mass-weighted center (V4)
    r_mass_center_mc = sqrt.((xxx .- xc_m).^2 .+ (yyy .- yc_m).^2 .+ (zzz .- zc_m).^2) |> collect
    beta_star_mc = compute_beta_star(r_mass_center_mc, rho,
        beta_star_r_min_kpc, beta_star_r_max_kpc, length_astro)

    # β* in tighter interval (robust to tidal-stripping artifact)
    beta_star_rm_inner = compute_beta_star(r_mass_center, rho,
        0.3, 0.5, length_astro)

    return DiagnosticSnapshot(
        t_current_code,
        rho_c, xc_m, yc_m, zc_m,
        bound_mass, interior_mass,
        beta_star_rm, beta_star_mc, beta_star_rm_inner,
        PE_abs, QE, KE,
    )
end

"""
$(TYPEDSIGNATURES)

Save the diagnostic time series to a JLD2 file and a CSV sidecar.

The file is written to
`joinpath(outputdir, "(title), (suffix) - Diagnostics.jld2")`
and contains the diagnostic DataFrame plus the per-shell velocity
field vectors at every saved timestep.
"""
function save_diagnostic_timeseries(
    diagnostics::Vector{DiagnosticSnapshot}, velocity_shells, outputdir, title, suffix)
    df = DataFrame(
        t          = [d.t                       for d in diagnostics],
        rho_c      = [d.rho_c                   for d in diagnostics],
        xc_mass    = [d.xc_mass                 for d in diagnostics],
        yc_mass    = [d.yc_mass                 for d in diagnostics],
        zc_mass    = [d.zc_mass                 for d in diagnostics],
        bound_mass = [d.bound_mass              for d in diagnostics],
        interior_mass = [d.interior_mass        for d in diagnostics],
        beta_star_rm         = [d.beta_star_rhomax         for d in diagnostics],
        beta_star_masscenter = [d.beta_star_masscenter for d in diagnostics],
        beta_star_rm_inner   = [d.beta_star_rhomax_interior for d in diagnostics],
        PE_abs     = [d.virial_V                 for d in diagnostics],
        QE         = [d.virial_Q                 for d in diagnostics],
        KE         = [d.virial_K                 for d in diagnostics],
    )
    save(joinpath(outputdir, "$(title), $(suffix) - Diagnostics.jld2"),
        Dict("df" => df, "velocity_shells" => velocity_shells))
    CSV.write(joinpath(outputdir, "$(title), $(suffix) - Diagnostics.csv"), df)
    return df
end

"""
$(TYPEDSIGNATURES)

Detect cuspy windows in a β*(t) time series.

A cuspy window is a connected interval of `t` over which
`β*(t) ≥ threshold` for at least `min_duration_code` (in code units).

Returns a NamedTuple with fields
- `windows::Vector{NamedTuple}` with `(t_start, t_end, t_peak, beta_max, duration)` per window
- `t_cuspy_frac::Float64` = total time in cuspy state / total time
"""
function detect_cuspy_windows(t, betastar;
    threshold = 0.5,
    min_duration_code = 0.05,
)
    if length(t) < 2
        return (windows = NamedTuple[], t_cuspy_frac = 0.0)
    end

    above = betastar .>= threshold
    windows = NamedTuple[]
    in_window = false
    t_start = 0.0
    beta_max = -Inf
    idx_start = 0
    idx_peak = 0

    for k in 1:length(t)
        if above[k]
            if !in_window
                in_window = true
                t_start = t[k]
                beta_max = betastar[k]
                idx_start = k
                idx_peak = k
            else
                if betastar[k] > beta_max
                    beta_max = betastar[k]
                    idx_peak = k
                end
            end
            if k == length(t) || !above[k+1]
                t_end = t[k]
                duration = t_end - t_start
                if duration >= min_duration_code
                    push!(windows, (
                        t_start = t_start,
                        t_end = t_end,
                        t_peak = t[idx_peak],
                        beta_max = beta_max,
                        duration = duration,
                        idx_start = idx_start,
                        idx_peak = idx_peak,
                    ))
                end
                in_window = false
            end
        end
    end

    t_cuspy_total = sum(w.duration for w in windows)
    t_total = t[end] - t[1]
    t_cuspy_frac = t_total > 0 ? t_cuspy_total / t_total : 0.0

    return (windows = windows, t_cuspy_frac = t_cuspy_frac)
end
