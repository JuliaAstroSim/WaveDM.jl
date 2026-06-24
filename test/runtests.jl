using Test
using FFTW
using Unitful, UnitfulAstro

using AstroSimBase
using PhysicalParticles
using AstroNbodySim

using WaveDM

using CairoMakie
CairoMakie.activate!(visible = false)


@time @testset "Unit tests" begin
    include("unit_test.jl")
end

@time @testset "Simulation" begin
    include("integration_test.jl")
end

@time @testset "3D SPE Validation Tests" begin
    include("simulation_SPE3D_validation.jl")

    # Run quick validation tests
    results = run_all_3d_validation_tests(verbose=false, quick=true)

    # Check that tests produce reasonable results
    @test haskey(results, "single_soliton")
    @test haskey(results, "soliton_collision")
    @test haskey(results, "soliton_binary")
    @test haskey(results, "convergence_spatial")
    @test haskey(results, "convergence_temporal")
    @test haskey(results, "energy_conservation")

    # Check energy conservation (should be reasonable for quick tests)
    @test results["single_soliton"][:energy_drift] < 1.0  # Less than 100% drift
    @test results["energy_conservation"][:energy_drift] < 1.0
end

@testset "Parallel backend" begin
    # The new `ParallelBackend`/`DeviceConfig` adaptive parallel system
    # is exercised end-to-end.  The distributed Poisson test is skipped
    # when fewer than 2 workers are available.
    include("test_parallel.jl")
end
