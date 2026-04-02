# FFT operators for spectral method

"""
$(TYPEDSIGNATURES)

Setup FFT operators for spectral method.
"""
function setup_fft_operators(Xmax, Ymax, Zmax, Nx, Ny, Nz, dt)
    kx = Vector{Float64}(LinRange(-Nx/4/Xmax, Nx/4/Xmax-1/2/Xmax, Nx))
    ky = Vector{Float64}(LinRange(-Ny/4/Ymax, Ny/4/Ymax-1/2/Ymax, Ny))
    kz = Vector{Float64}(LinRange(-Nz/4/Zmax, Nz/4/Zmax-1/2/Zmax, Nz))
    Laplacian = [(2π*im*kx[i])^2 + (2π*im*ky[j])^2 + (2π*im*kz[k])^2 for i in 1:Nx, j in 1:Ny, k in 1:Nz]
    linear_phase = fftshift(exp.(im * Laplacian * dt / 2))
    return linear_phase
end

"""
$(TYPEDSIGNATURES)

Setup absorption boundary conditions.
"""
function setup_absorption_boundary(Xmax, Ymax, Zmax, x, y, z, absorb_coeff, dt)
    wx = Xmax/50
    wy = Ymax/50
    wz = Zmax/50
    return exp.(-absorb_coeff*(6 .- tanh.((x.+Xmax)./wx) .+ tanh.((x.-Xmax)./wx) .- tanh.((y.+Ymax)./wy) .+ tanh.((y.-Ymax)./wy) .- tanh.((z.+Zmax)./wz) .+ tanh.((z.-Zmax)./wz))*dt)
end
