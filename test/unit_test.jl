@testset "setup_grid" begin
    Xmax = Ymax = Zmax = 2.0
    Nx = Ny = Nz = 64
    
    x, y, z, Δ, unit_cell_volumn = WaveDM.setup_grid(Xmax, Ymax, Zmax, Nx, Ny, Nz)
    
    @test minimum(x) ≈ -Xmax
    @test maximum(x) ≈ Xmax
    @test length(x) == Nx
    @test length(y) == Ny
    @test length(z) == Nz
    
    @test length(Δ) == 3
    @test Δ[1] ≈ (2*Xmax)/(Nx-1)
    @test Δ[2] ≈ (2*Ymax)/(Ny-1)
    @test Δ[3] ≈ (2*Zmax)/(Nz-1)
    
    expected_vol = prod(Δ)
    @test unit_cell_volumn ≈ expected_vol
    
    x8, y8, z8, Δ8, _ = WaveDM.setup_grid(1.0, 1.0, 1.0, 8, 8, 8)
    @test length(x8) == 8
    @test Δ8[1] ≈ (2*1.0)/(8-1)
end

@testset "setup_coordinates" begin
    Xmax = Ymax = Zmax = 2.0
    Nx = Ny = Nz = 9
    x = collect(LinRange(-Xmax, Xmax, Nx))
    y = collect(LinRange(-Ymax, Ymax, Ny))
    z = collect(LinRange(-Zmax, Zmax, Nz))
    oneMatrix = ones(Nx, Ny, Nz)
    
    xxx, yyy, zzz, r = WaveDM.setup_coordinates(x, y, z, Nx, Ny, Nz, oneMatrix; DA = collect)
    
    @test size(xxx) == (Nx, Ny, Nz)
    @test size(yyy) == (Nx, Ny, Nz)
    @test size(zzz) == (Nx, Ny, Nz)
    @test size(r) == (Nx, Ny, Nz)
    
    center_idx = div(Nx, 2) + 1
    @test xxx[center_idx, center_idx, center_idx] ≈ 0.0
    @test yyy[center_idx, center_idx, center_idx] ≈ 0.0
    @test zzz[center_idx, center_idx, center_idx] ≈ 0.0
    
    @test r[center_idx, center_idx, center_idx] ≈ 0.0 atol = 1e-10
    @test all(r .>= 0)
    
    @test xxx[:, 1, 1] == x
    @test yyy[1, :, 1] == y
    @test zzz[1, 1, :] == z
end

@testset "solve_vector_equation" begin
    Nx = Ny = Nz = 8
    vx = rand(Nx, Ny, Nz)
    vy = rand(Nx, Ny, Nz)
    vz = rand(Nx, Ny, Nz)
    
    Δ = [0.25, 0.25, 0.25]
    
    θ = WaveDM.solve_vector_equation(collect(vx), collect(vy), collect(vz), Δ...)
    
    @test length(θ) == Nx * Ny * Nz
    @test !any(isnan.(θ))
    @test !any(isinf.(θ))
    
    phase = exp.(im*θ)
    @test all(abs.(phase) .≈ 1.0)
end
