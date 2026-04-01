# SPE schrodinger module

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
$(TYPEDSIGNATURES)

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
