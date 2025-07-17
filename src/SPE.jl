using PhysicalParticles.NumericalIntegration

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

# function solve_vector_equation(vx, vy, vz, Δx, Δy, Δz)
#     Nx, Ny, Nz = size(vx)
#     θ = zeros(Nx, Ny, Nz)
#     for i in 2:Nx
#         θ[i,:,:] = θ[i-1,:,:] + Δx * vx[i-1,:,:]
#     end
#     for j in 2:Ny
#         θ[:,j,:] = θ[:,j-1,:] + Δy * vy[:,j-1,:]
#     end
#     for k in 2:Nz
#         θ[:,:,k] .= θ[:,:,k-1] .+ Δz * vz[:,:,k-1]
#     end
#     return θ
# end