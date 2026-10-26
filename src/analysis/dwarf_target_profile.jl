# Dwarf-galaxy target halo profile analysis
#
# Extracted from the SPE-WaveDM example scripts:
#   * `SPE_WaveDM_dwarfs_paper_massive.jl`
#   * `SPE_WaveDM_dwarfs_paper_massive_compat.jl`


"""
$(TYPEDSIGNATURES)

Extract the UFD (`Zhao`-profile) target parameters for galaxy index `i`
from a `df_UFDs::DataFrame` produced by `AstroIC.load_data_UFDs`.

Returns a `NamedTuple` with the base value and upper/lower perturbations
for every Zhao parameter.
"""
function get_dwarf_UFDs_params(i::Integer, df_UFDs)
    return (
        ρ0   = df_UFDs.rho0[i],
        ρ0_u = df_UFDs.rho0_u[i],
        ρ0_d = df_UFDs.rho0_d[i],
        rs   = df_UFDs.b[i],
        rs_u = df_UFDs.b_u[i],
        rs_d = df_UFDs.b_d[i],
        α    = df_UFDs.alpha[i],
        α_u  = df_UFDs.alpha[i] + df_UFDs.alpha_u[i],
        α_d  = df_UFDs.alpha[i] + df_UFDs.alpha_d[i],
        β    = df_UFDs.beta[i],
        β_u  = df_UFDs.beta[i] + df_UFDs.beta_u[i],
        β_d  = df_UFDs.beta[i] + df_UFDs.beta_d[i],
        γ    = df_UFDs.gamma[i],
        γ_u  = df_UFDs.gamma[i] + df_UFDs.gamma_u[i],
        γ_d  = df_UFDs.gamma[i] + df_UFDs.gamma_d[i],
    )
end

"""
$(TYPEDSIGNATURES)

Extract the gNFW target parameters for galaxy index `i` from a
`dfDwarf::DataFrame`.  Returns base value and upper/lower perturbations
for `(β, ρ0, rs)` and the virial-region bounds (in kpc) used for the
inner profile-slope fit.
"""
function get_dwarf_gNFW_params(i::Integer, dfDwarf)
    return (
        β    = dfDwarf.beta[i],
        β_u  = dfDwarf.beta[i] + dfDwarf.beta_u[i],
        β_d  = dfDwarf.beta[i] - dfDwarf.beta_d[i],
        ρ0   = dfDwarf.rho0[i],
        ρ0_u = dfDwarf.rho0[i] + dfDwarf.rho0_u[i],
        ρ0_d = dfDwarf.rho0[i] - dfDwarf.rho0_d[i],
        rs   = dfDwarf.rs[i],
        rs_u = dfDwarf.rs[i] + dfDwarf.rs_u[i],
        rs_d = dfDwarf.rs[i] - dfDwarf.rs_d[i],
        inner_r_vir_min = ustrip(u"kpc", dfDwarf.R_vir1[i]),
        inner_r_vir_max = ustrip(u"kpc", dfDwarf.R_vir2[i]),
    )
end

"""
$(TYPEDSIGNATURES)

Extract NFW target parameters for a X-ray ETG (only `(ρ0, rs)` are needed
for that model, no error bands).
"""
function get_Xray_ETG_NFW_params(i::Integer, dfSPARC_Xray_ETGs)
    return (
        ρ0 = dfSPARC_Xray_ETGs.rho0[i],
        rs = dfSPARC_Xray_ETGs.rs[i],
    )
end

"""
$(TYPEDSIGNATURES)

Dispatcher that returns a `NamedTuple` of the target-profile parameters
for galaxy index `Galaxy_i` and `target_profile_model` ∈
{`:dwarf_Zhao`, `:dwarf_UFDs`, `:dwarf_gNFW`, `:dwarf_massive`,
`:dwarf_NFW`, `:Xray_ETG_NFW`}.

The relevant data tables must be passed as keyword arguments:
- `df_UFDs`            – required for `:dwarf_Zhao` / `:dwarf_UFDs`
- `dfDwarf`            – required for `:dwarf_gNFW` / `:dwarf_massive`
- `dfSPARC_Xray_ETGs`  – required for `:Xray_ETG_NFW`

Any other symbol raises an `ArgumentError`.
"""
function get_target_profile_params(Galaxy_i::Integer, target_profile_model::Symbol;
                                   df_UFDs           = nothing,
                                   dfDwarf           = nothing,
                                   dfSPARC_Xray_ETGs = nothing)
    if target_profile_model == :dwarf_Zhao || target_profile_model == :dwarf_UFDs
        df_UFDs === nothing && error("`:dwarf_Zhao`/`:dwarf_UFDs` requires `df_UFDs`")
        return get_dwarf_UFDs_params(Galaxy_i, df_UFDs)
    elseif target_profile_model == :dwarf_gNFW || target_profile_model == :dwarf_massive
        dfDwarf === nothing && error("`:dwarf_gNFW`/`:dwarf_massive` requires `dfDwarf`")
        return get_dwarf_gNFW_params(Galaxy_i, dfDwarf)
    elseif target_profile_model == :dwarf_NFW || target_profile_model == :Xray_ETG_NFW
        dfSPARC_Xray_ETGs === nothing &&
            error("`:dwarf_NFW`/`:Xray_ETG_NFW` requires `dfSPARC_Xray_ETGs`")
        return get_Xray_ETG_NFW_params(Galaxy_i, dfSPARC_Xray_ETGs)
    else
        throw(ArgumentError("Unsupported target_profile_model :$(target_profile_model)"))
    end
end

# ------------------------------------------------------------------------------
# Target profile evaluation + error bands
# ------------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Evaluate a target halo profile on `sample_r` (in `kpc`) and return the
density in code units (the unit of `density_astro`).

The function dispatches on the value of `model`.  For Zhao/gNFW models
the same `params` NamedTuple returned by `get_target_profile_params`
is expected.
"""
function evaluate_target_profile(model::Symbol, params, sample_r, density_astro)
    sample_r_q = sample_r * u"kpc"
    if model == :dwarf_Zhao || model == :dwarf_UFDs
        m = Zhao(params.ρ0, params.rs, params.α, params.β, params.γ)
        ρ = GalacticDynamics.density.(m, sample_r_q) ./ density_astro
    elseif model == :dwarf_gNFW || model == :dwarf_massive
        m = gNFW(params.β, params.ρ0, params.rs)
        ρ = GalacticDynamics.density.(m, sample_r_q) ./ density_astro
    elseif model == :dwarf_NFW || model == :Xray_ETG_NFW
        m = NFW(params.ρ0, params.rs)
        ρ = GalacticDynamics.density.(m, sample_r_q) ./ density_astro
    else
        throw(ArgumentError("Unsupported target_profile_model :$(model)"))
    end
    return upreferred.(ρ)
end

"""
$(TYPEDSIGNATURES)

Compute the upper / lower error bands of a target halo density profile
by taking the Cartesian product of `(base, +u, +d)` perturbations of
each parameter and computing the per-radius `extrema` over all
``3^n`` models.  This matches the logic that was inlined in the
SPE-WaveDM dwarf-massive scripts.

# Returns a `NamedTuple`
- `ρ_halo_target_u` / `ρ_halo_target_d`: vectors of upper/lower envelope
  densities in code units (same length as `sample_r`).
- `densities`: a `(length(sample_r), length(models))` matrix of all
  evaluated densities (in `Msun/kpc^3`) — useful for additional plotting.
- `models`: the list of constructed `GalacticDynamics` model objects.
"""
function compute_target_density_with_errors(model::Symbol, params, sample_r, density_astro)
    sample_r_q = sample_r * u"kpc"

    if model == :dwarf_Zhao || model == :dwarf_UFDs
        param_sets = Base.Iterators.product(
            [params.ρ0, params.ρ0_u, params.ρ0_d],
            [params.rs, params.rs_u, params.rs_d],
            [params.α,  params.α_u,  params.α_d],
            [params.β,  params.β_u,  params.β_d],
            [params.γ,  params.γ_u,  params.γ_d],
        )
        models = [Zhao(p...) for p in param_sets]
    elseif model == :dwarf_gNFW || model == :dwarf_massive
        param_sets = Base.Iterators.product(
            [params.β,  params.β_u,  params.β_d],
            [params.ρ0, params.ρ0_u, params.ρ0_d],
            [params.rs, params.rs_u, params.rs_d],
        )
        models = [gNFW(p...) for p in param_sets]
    elseif model == :dwarf_NFW || model == :Xray_ETG_NFW
        # No error bands available for the bare NFW X-ray ETG profile.
        param_sets = Base.Iterators.product([params.ρ0], [params.rs])
        models = [NFW(p...) for p in param_sets]
    else
        throw(ArgumentError("Unsupported target_profile_model :$(model)"))
    end

    densities = zeros(length(sample_r), length(models))
    ρ_halo_target_u = similar(sample_r, Float64)
    ρ_halo_target_d = similar(sample_r, Float64)

    for i in eachindex(sample_r)
        densities[i, :] .= ustrip.(u"Msun/kpc^3",
            GalacticDynamics.density.(models, sample_r_q[i]))
        u, d = extrema(densities[i, :])
        ρ_halo_target_u[i] = u
        ρ_halo_target_d[i] = d
    end

    return (;
        ρ_halo_target_u   = ρ_halo_target_u ./ ustrip(u"Msun/kpc^3", density_astro),
        ρ_halo_target_d   = ρ_halo_target_d ./ ustrip(u"Msun/kpc^3", density_astro),
        densities         = densities,
        models            = models,
    )
end

# ------------------------------------------------------------------------------
# Virial radius
# ------------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Return the radius (in `kpc`) at which `ρ_halo_target` (in code units)
crosses the cosmic critical density `3 H² / (8 π G)`, or `nothing` if
no crossing is found.

`sample_r` is in `kpc`; `density_astro` is the WaveDM density unit used
to convert `ρ_halo_target` into physical `Msun/kpc^3` units.
"""
function find_r_virial(ρ_halo_target, sample_r, density_astro;
                       H0::Real=C.H, G::Real=C.G)
    rho_critical_phys = uconvert(u"Msun/kpc^3", 3 * H0^2 / (8π * G))
    rho_critical_code = upreferred(rho_critical_phys / density_astro)
    id_critical = findfirstvalue(ρ_halo_target, rho_critical_code)
    return id_critical === nothing ? nothing : sample_r[id_critical]
end


"""
$(TYPEDSIGNATURES)

Lelli 2017 (Eq. 11) empirical RAR:

    g_obs = g_bar / (1 - exp(-sqrt(g_bar / g†)))

with the only free parameter being `g†`.  At `g_bar ≫ g†` this
asymptotes to `g_obs ≈ g_bar` (Newtonian / baryon-dominated), and at
`g_bar ≪ g†` it tends to `g_obs ≈ sqrt(g_bar · g†)` (MOND
asymptote).  The fiducial fit gives `g† = 1.20 × 10⁻¹⁰ m s⁻²`.

Accepts either bare-Float64 (assumed `m/s²`) or `Unitful` quantities.
"""
function Lelli2017_RAR(gbar, gdagger::Real = 1.20e-10)
    gbar_q = if gbar isa Quantity
        ustrip(u"m/s^2", gbar)
    else
        float(gbar)
    end
    return gbar_q / (1.0 - exp(-sqrt(gbar_q / gdagger)))
end

"""
$(TYPEDSIGNATURES)

Inverse of [`Lelli2017_RAR`](@ref) — given `g_obs`, solve for `g_bar`
with `gdagger` fixed.  Returns a `Float64` in `m/s²`.  For
`g_obs < gdagger` the mapping is two-valued; we return the
physically relevant (lower) root `g_bar < g_obs`.

The Lelli 2017 Eq. (11) functional form has poor linear conditioning
(g_obs ≈ g_bar in the Newtonian limit, g_obs ≈ sqrt(g_bar · gdagger)
in the MOND limit), so we solve

    log10(RAR(g_bar)) = log10(g_obs)

in log-space using Newton's method with the exact derivative

    d log RAR / d log g_bar = (g_bar / RAR) · dRAR/dg_bar

`dRAR/dg_bar = 1/(1 - e^{-s}) - g_bar · e^{-s} / (2 sqrt(g_bar·gdagger) · (1 - e^{-s})²)`

where `s = sqrt(g_bar / gdagger)`.  Convergence to `|Δ log10(g_bar)| < 1e-10`
takes 2-3 iterations for any `g_obs ∈ [1e-15, 1e-5] m/s²`.
"""
function Lelli2017_RAR_inv(gobs, gdagger::Real = 1.20e-10)
    gobs_q = if gobs isa Quantity
        ustrip(u"m/s^2", gobs)
    else
        float(gobs)
    end
    gobs_q <= 0 && return 0.0

    # f(u) = log10(RAR(10^u)) - log10(g_obs),  f'(u) = d(log RAR)/d(log g_bar)
    function f_df(u)
        gb     = 10.0^u
        s      = sqrt(gb / gdagger)
        h      = 1.0 - exp(-s)
        rar    = gb / h
        f_     = log10(rar) - log10(gobs_q)
        dsdb   = 1.0 / (2.0 * sqrt(gb * gdagger))
        drar   = 1.0 / h - gb * exp(-s) * dsdb / h^2
        df_    = (gb / rar) * drar
        return f_, df_
    end

    # Initial guess: invert RAR(g_bar) ≈ g_bar / (1 - exp(-sqrt(g_bar/gdagger)))
    # ≈ g_obs exactly gives log10(g_bar) = log10(g_obs) + log10(1 - exp(-sqrt(g_bar/gdagger))).
    # Use the g_obs ≈ sqrt(g_bar · gdagger) MOND limit as the seed.
    s0 = sqrt(max(gobs_q, 1e-30) / gdagger)
    u0 = log10(gobs_q) - log10(1.0 - exp(-s0))

    # Newton iterations (capped at 20 — converges in 2-3 for any
    # well-formed g_obs).
    for _ in 1:20
        f1, df1 = f_df(u0)
        if abs(f1) < 1e-12 || abs(df1) < 1e-30
            break
        end
        u0 = u0 - f1 / df1
    end
    return 10.0^u0
end
