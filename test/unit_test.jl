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

@testset "Lelli 2017 RAR" begin
    using WaveDM: Lelli2017_RAR, Lelli2017_RAR_inv, a0_Lelli2017
    using AstroIC: lelli2017_gdagger_mps2

    # The published fiducial value
    @test a0_Lelli2017 == 1.20e-10
    @test lelli2017_gdagger_mps2() == 1.20e-10 * u"m/s^2"

    # At gbar >> g†, the RAR asymptotes to gbar (Newtonian limit)
    gbar_high = 1e-8                       # 100 × g†
    @test Lelli2017_RAR(gbar_high) ≈ gbar_high  rtol = 1e-3

    # At gbar == g†, RAR ≈ g† / (1 - exp(-1)) ≈ 1.582 × g†
    gbar_eq  = a0_Lelli2017
    g_obs_eq = Lelli2017_RAR(gbar_eq)
    @test g_obs_eq ≈ a0_Lelli2017 / (1 - exp(-1))  rtol = 1e-6

    # At gbar << g†, RAR ≈ sqrt(gbar · g†) (MOND asymptote)
    gbar_low = 1e-13
    @test Lelli2017_RAR(gbar_low) ≈ sqrt(gbar_low * a0_Lelli2017)  rtol = 5e-3

    # Round-trip: g_obs → gbar → g_obs should reproduce within 1e-6
    for gbar0 in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8)
        gobs = Lelli2017_RAR(gbar0)
        gbar_back = Lelli2017_RAR_inv(gobs)
        @test gbar_back ≈ gbar0  rtol = 1e-6
    end

    # Unitful input
    @test Lelli2017_RAR(1e-10 * u"m/s^2") ≈ Lelli2017_RAR(1e-10)
    @test Lelli2017_RAR_inv(1e-10 * u"m/s^2") ≈ Lelli2017_RAR_inv(1e-10)
end

@testset "build_xray_etg_halo (flat branch, analytic)" begin
    # Re-implement the relevant calls directly (no WaveDM dependency on
    # the example script's helpers).  Test the analytic flat profile
    # branch with hand-crafted values.
    using WaveDM: build_flat_halo_density
    r_kpc = collect(LinRange(0.1, 30, 1000))
    M_total = 1e12
    r_flat = 5.0
    rho_min = 1e6
    R_max = 30.0
    rho_kpc3, rho_flat = build_flat_halo_density(r_kpc, M_total, r_flat, rho_min, R_max)
    @test length(rho_kpc3) == length(r_kpc)
    @test all(rho_kpc3 .>= 0)
    @test all(rho_kpc3[findall(r_kpc .<= r_flat)] .≈ rho_flat)
    @test all(rho_kpc3[findall(r_flat .< r_kpc .<= R_max)] .≈ rho_min)
    @test all(rho_kpc3[findall(r_kpc .> R_max)] .== 0)

    # Invalid: R_max must exceed r_flat
    @test_throws ErrorException build_flat_halo_density(r_kpc, M_total, r_flat, rho_min, r_flat)
    @test_throws ErrorException build_flat_halo_density(r_kpc, M_total, -1.0, rho_min, R_max)
end
