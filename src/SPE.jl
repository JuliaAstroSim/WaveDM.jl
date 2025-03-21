function solve_vector_equation(vx, vy, vz, Δx, Δy, Δz)
    Nx, Ny, Nz = size(vx)
    θ = zeros(Nx, Ny, Nz)
    for i in 2:Nx
        θ[i,:,:] = θ[i-1,:,:] + Δx * vx[i-1,:,:]
    end
    for j in 2:Ny
        θ[:,j,:] = θ[:,j-1,:] + Δy * vy[:,j-1,:]
    end
    for k in 2:Nz
        θ[:,:,k] .= θ[:,:,k-1] .+ Δz * vz[:,:,k-1]
    end
    return θ
end