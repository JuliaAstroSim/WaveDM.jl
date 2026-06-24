"""
$(TYPEDSIGNATURES)

Compute the cartesian probability-current velocity field
`(vx, vy, vz)` and the canonical-velocity field
`(wx, wy, wz)` of a complex wave function `ψ`.

# Arguments
- `ψ::AbstractArray{<:Complex}`: 3-D complex wave function on a uniform
  grid of spacing `Δ = (Δx, Δy, Δz)`.
- `Δ`: tuple `(Δx, Δy, Δz)` of grid spacings.
- `ρ = abs2.(ψ)`: pre-computed density.  Passing it avoids recomputing
  `abs2.(ψ)` and lets the caller reuse the value.

# Returns
NamedTuple with `vx, vy, vz, wx, wy, wz` — each a real `Array` of the
same shape as `ψ`.  NaNs in `ρ` (e.g. at the very centre of the grid)
propagate into the returned velocities.
"""
function compute_waveDM_velocity_field(ψ::AbstractArray{<:Complex},
                                       Δ,
                                       ρ::AbstractArray = abs2.(ψ))
    ψx, ψy, ψz   = grad_central(Δ..., ψ)
    ψc           = conj.(ψ)
    ψcx, ψcy, ψcz = grad_central(Δ..., ψc)

    ρvx = real.((ψc .* ψx .- ψ .* ψcx) ./ (2im))
    ρwx = real.((ψc .* ψx .+ ψ .* ψcx) ./ 2)

    ρvy = real.((ψc .* ψy .- ψ .* ψcy) ./ (2im))
    ρwy = real.((ψc .* ψy .+ ψ .* ψcy) ./ 2)

    ρvz = real.((ψc .* ψz .- ψ .* ψcz) ./ (2im))
    ρwz = real.((ψc .* ψz .+ ψ .* ψcz) ./ 2)

    return (
        vx = ρvx ./ ρ,
        vy = ρvy ./ ρ,
        vz = ρvz ./ ρ,
        wx = ρwx ./ ρ,
        wy = ρwy ./ ρ,
        wz = ρwz ./ ρ,
    )
end

"""
$(TYPEDSIGNATURES)

Rotate a cartesian velocity field to spherical components around
`rho_max_id`.  Inputs are the displacements of every grid point from
the density maximum.

By default, returns the `v` field only; pass `with_w=true` to also
return the canonical-velocity `w` field rotated to spherical.
"""
function cartesian_to_spherical_velocity_field(vx, vy, vz, wx, wy, wz,
                                                xxx, yyy, zzz, r_mass_center,
                                                rho_max_id;
                                                with_w::Bool=true)
    x_shift = xxx .- xxx[rho_max_id]
    y_shift = yyy .- yyy[rho_max_id]
    z_shift = zzz .- zzz[rho_max_id]

    vr, vθ, vϕ = vec_cartesian_to_spherical(vx, vy, vz,
                                            x_shift, y_shift, z_shift,
                                            r_mass_center)
    if with_w
        wr, wθ, wϕ = vec_cartesian_to_spherical(wx, wy, wz,
                                                x_shift, y_shift, z_shift,
                                                r_mass_center)
        return (; vr, vθ, vϕ, wr, wθ, wϕ)
    else
        return (; vr, vθ, vϕ)
    end
end

"""
$(TYPEDSIGNATURES)

Bin the spherical velocity components `(vr, vθ, vϕ)` on a radial grid
`r_mass_center` and return their standard deviations `σr, σθ, σϕ` in
each spherical shell.
"""
function compute_velocity_dispersion_profile(r_mass_center, vr, vθ, vϕ;
                                             r_filter=nothing,
                                             section::Integer=100,
                                             uniform_interval::Bool=true)
    if r_filter === nothing
        r_filter = trues(length(r_mass_center))
    end

    r_mean, vr_mean, r_std, σr = distribution(r_mass_center[r_filter], vr[r_filter];
                                              section, uniform_interval)
    _,      vθ_mean, _,     σθ = distribution(r_mass_center[r_filter], vθ[r_filter];
                                              section, uniform_interval)
    _,      vϕ_mean, _,     σϕ = distribution(r_mass_center[r_filter], vϕ[r_filter];
                                              section, uniform_interval)

    return (
        r_mean = r_mean,
        vr_mean = vr_mean,
        vθ_mean = vθ_mean,
        vϕ_mean = vϕ_mean,
        σr = σr,
        σθ = σθ,
        σϕ = σϕ,
    )
end
