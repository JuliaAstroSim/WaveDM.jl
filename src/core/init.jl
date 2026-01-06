# SPE initialization module

"""
$(TYPEDSIGNATURES)

Solve vector equation ∇θ = v for phase θ given velocity field v.
"""
function solve_vector_equation(vx, vy, vz, Δx, Δy, Δz)
    Nx, Ny, Nz = size(vx)
    θ = zeros(Nx, Ny, Nz)
    x = collect(0:Nx-1) * Δx
    y = collect(0:Ny-1) * Δy
    z = collect(0:Nz-1) * Δz
    for k in 1:Nz
        for j in 1:Ny
            θ[:,j,k] .+= cumul_integrate(x, vx[:,j,k])
        end
    end
    for k in 1:Nz
        for i in 1:Nx
            θ[i,:,k] .+= cumul_integrate(y, vy[i,:,k])
        end
    end
    for j in 1:Ny
        for i in 1:Nx
            θ[i,j,:] .+= cumul_integrate(z, vz[i,j,:])
        end
    end
    return θ
end

"""
    setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)

Setup computational grid and coordinates.
"""
function setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)
    x = collect(LinRange(-Xmax, Xmax, Nx))
    y = collect(LinRange(-Ymax, Ymax, Ny))
    z = collect(LinRange(-Zmax, Zmax, Nz))
    Δ = [x[2]-x[1], y[2]-y[1], z[2]-z[1]]
    unit_cell_volumn = prod(Δ)
    return x, y, z, Δ, unit_cell_volumn
end

"""
    setup_coordinates(x, y, z, Nx, Ny, Nz; DA = collect)

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
    setup_initial_conditions(IC, xxx, yyy, zzz; DA = collect)

Setup initial wave function.
"""
function setup_initial_conditions(IC, xxx, yyy, zzz; DA = collect)
    if IC isa Function
        ψ = ComplexF64.(IC.(xxx, yyy, zzz))
    else
        ψ = DA(IC)
    end
    sqrt_rho = abs.(ψ)
    rho = sqrt_rho.^2
    return ψ, sqrt_rho, rho
end

"""
    compute_timestep(Δx, Φ_all, κ, ψ, Tmax, Nt, autoset_timestep, autoset_timestep_ratio)

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
    
    t = collect(LinRange(0, Tmax, Nt))
    dt = Float32(Tmax / Nt)
    
    return t, dt, Nt
end

"""
    setup_fft_operators(Xmax, Ymax, Zmax, Nx, Ny, Nz, dt)

Setup FFT operators for spectral method.
"""
function setup_fft_operators(Xmax, Ymax, Zmax, Nx, Ny, Nz, dt)
    kx = collect(LinRange(-Nx/4/Xmax, Nx/4/Xmax-1/2/Xmax, Nx))
    ky = collect(LinRange(-Ny/4/Ymax, Ny/4/Ymax-1/2/Ymax, Ny))
    kz = collect(LinRange(-Nz/4/Zmax, Nz/4/Zmax-1/2/Zmax, Nz))
    Laplacian = [(2π*im*kx[i])^2 + (2π*im*ky[j])^2 + (2π*im*kz[k])^2 for i in 1:Nx, j in 1:Ny, k in 1:Nz]
    linear_phase = fftshift(exp.(im * Laplacian * dt / 2))
    return linear_phase
end

"""
    setup_absorption_boundary(Xmax, Ymax, Zmax, x, y, z, absorb_coeff, dt)

Setup absorption boundary conditions.
"""
function setup_absorption_boundary(Xmax, Ymax, Zmax, x, y, z, absorb_coeff, dt)
    wx = Xmax/50
    wy = Ymax/50
    wz = Zmax/50
    return exp.(-absorb_coeff*(6 .- tanh.((x.+Xmax)./wx) .+ tanh.((x.-Xmax)./wx) .- tanh.((y.+Ymax)./wy) .+ tanh.((y.-Ymax)./wy) .- tanh.((z.+Zmax)./wz) .+ tanh.((z.-Zmax)./wz))*dt)
end
