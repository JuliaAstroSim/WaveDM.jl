# SPE poisson module

"""
$(TYPEDSIGNATURES)

Setup computational grid and coordinates.
"""
function setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)
    x = Vector{Float64}(LinRange(-Xmax, Xmax, Nx))
    y = Vector{Float64}(LinRange(-Ymax, Ymax, Ny))
    z = Vector{Float64}(LinRange(-Zmax, Zmax, Nz))
    Δ = [x[2]-x[1], y[2]-y[1], z[2]-z[1]]
    unit_cell_volumn = prod(Δ)
    return x, y, z, Δ, unit_cell_volumn
end

"""
$(TYPEDSIGNATURES)

Setup coordinate arrays.
"""
function setup_coordinates(x, y, z, Nx, Ny, Nz, oneMatrix = ones(Nx, Ny, Nz); DA = collect)
    xxx = oneMatrix .* x |> DA
    yyy = oneMatrix .* y' |> DA
    zzz = oneMatrix .* reshape(z, 1, 1, Nz) |> DA
    r = sqrt.(xxx.^2 + yyy.^2 + zzz.^2)
    return xxx, yyy, zzz, r
end

"""
$(TYPEDSIGNATURES)

Compute appropriate timestep.
"""
function compute_timestep(Δx, Φ_all, κ, ψ, Tmax, Nt, autoset_timestep, autoset_timestep_ratio)
    Δt_phase = 2 * Δx^2 / π
    Δt_grav = 2π / maximum(abs.(Φ_all))
    Δt_SI = 2π / abs(κ) / maximum(abs.(ψ))
    @info "Suggested Δt less than [2(Δx)²/π, 2π/max{Φ}, 2π/(|κ| max{ψ})] = [$(Δt_phase), $(Δt_grav), $(Δt_SI)]"

    if autoset_timestep
        Δt_limit = min(Δt_phase, Δt_grav, Δt_SI)
        Nt = ceil(Int, Tmax / (autoset_timestep_ratio * Δt_limit))
    end

    t = Vector{Float64}(LinRange(0, Tmax, Nt))
    dt = Float32(Tmax / Nt)

    return t, dt, Nt
end