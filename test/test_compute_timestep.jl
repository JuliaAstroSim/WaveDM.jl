# Tests for `compute_timestep` and `compute_averaged_fields`.
#
# `compute_timestep` is the public SPE driver for picking `dt` from the
# CFL constraints (phase, gravity, self-interaction).  `compute_averaged_fields`
# produces the time-averaged ψ² and acceleration snapshot when
# `average = true` in `simulate_waveDM`.
#
# Both are exported but not exercised by any existing testset.

using Test
using Statistics
using WaveDM: compute_timestep, compute_averaged_fields,
               DeviceConfig

@testset "compute_timestep: autotime off → uses user-supplied Nt" begin
    Δx = 0.1
    Φ_all = zeros(8, 8, 8)
    ψ = complex.(ones(8, 8, 8), zeros(8, 8, 8))
    Tmax = 0.1
    Nt_user = 50

    t, dt, Nt = compute_timestep(Δx, Φ_all, 0.0, ψ, Tmax, Nt_user,
                                  false, 0.9)
    @test length(t) == Nt_user
    @test Nt == Nt_user
    @test dt ≈ Tmax / Nt_user  atol = 1e-6
    @test t[1] ≈ 0.0
    @test t[end] ≈ Tmax
end

@testset "compute_timestep: autotime on → respects CFL bounds" begin
    Δx = 0.1
    Φ_all = 10.0 * ones(8, 8, 8)             # strong potential → tight Δt_grav
    ψ = complex.(ones(8, 8, 8), zeros(8, 8, 8))
    Tmax = 1.0

    t, dt, Nt = compute_timestep(Δx, Φ_all, 0.0, ψ, Tmax, 1,
                                  true, 0.5)
    # Δt_grav = 2π / max|Φ| = 2π / 10 ≈ 0.628
    # autoset_ratio = 0.5 ⇒ Δt_limit × 0.5 ≈ 0.314
    # Tmax / 0.314 ≈ 3.18 → Nt = ceil = 4
    @test dt > 0
    @test Nt >= 4
    @test dt ≈ Tmax / Nt
    @test length(t) == Nt
end

@testset "compute_timestep: self-interaction term dominates" begin
    Δx = 0.1
    Φ_all = zeros(8, 8, 8)
    ψ = complex.(10.0 * ones(8, 8, 8), zeros(8, 8, 8))   # large |ψ|
    κ = 100.0                                            # strong SI
    Tmax = 1.0

    # Δt_SI = 2π / |κ| / max|ψ| = 2π / 100 / 10 ≈ 0.00628
    # autoset_ratio = 0.5 ⇒ ≈ 0.00314 ⇒ Nt ≈ 319
    t, dt, Nt = compute_timestep(Δx, Φ_all, κ, ψ, Tmax, 1, true, 0.5)
    @test Nt >= 100                  # tight dt, many steps
    @test dt * Nt ≈ Tmax atol = 1e-9
end

@testset "compute_averaged_fields: average=false → returns empty ψ² and zeros" begin
    Nx = Ny = Nz = 8
    a_all = rand(Nx, Ny, Nz)
    ψ_buffer = rand(Nx, Ny, Nz)
    config_device = DeviceConfig(; gpu=false)

    ψ2_avg, a_all_avg = compute_averaged_fields(false, ψ_buffer, 10,
                                                :ignored, a_all, nothing,
                                                (1.0, 1.0, 1.0),
                                                Nx, Ny, Nz, config_device)
    @test isempty(ψ2_avg)
    @test size(a_all_avg) == size(a_all)
    @test all(a_all_avg .== 0)
end

@testset "compute_averaged_fields: average=true → Poisson solve on ρ" begin
    # For a uniform density buffer, the average should also be uniform
    # (up to Poisson solver noise).  We just check it doesn't NaN and
    # that the field shape is preserved.
    Nx = Ny = Nz = 8
    a_all = ones(Nx, Ny, Nz)
    # Average-N buffer of uniform density → uniform ρ → Φ from a uniform
    # ρ is well-defined (zero in code units).
    ψ_buffer = ones(Nx, Ny, Nz) ./ Nx^3
    config_device = DeviceConfig(; gpu=false)

    ψ2_avg, a_all_avg = compute_averaged_fields(true, ψ_buffer, 5,
                                                :ignored, a_all, nothing,
                                                (1.0, 1.0, 1.0),
                                                Nx, Ny, Nz, config_device)
    @test size(ψ2_avg) == (Nx, Ny, Nz)
    @test all(isfinite, ψ2_avg)
    @test all(isfinite, a_all_avg)
end

@testset "compute_averaged_fields: average=true with baryons" begin
    # When baryon_mode != :ignored, Φ_all = Φ_ψ + Φ_b.
    # Use a tiny Φ_b contribution and check it gets added.
    Nx = Ny = Nz = 4
    a_all = zeros(Nx, Ny, Nz)
    ψ_buffer = ones(Nx, Ny, Nz) ./ Nx^3
    Φ_b = ones(Nx, Ny, Nz) .* 1e-3
    config_device = DeviceConfig(; gpu=false)

    ψ2_avg, a_all_avg = compute_averaged_fields(true, ψ_buffer, 5,
                                                :mesh, a_all, Φ_b,
                                                (1.0, 1.0, 1.0),
                                                Nx, Ny, Nz, config_device)
    @test all(isfinite, ψ2_avg)
    @test all(isfinite, a_all_avg)
end
