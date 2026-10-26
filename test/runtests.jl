using Test
using FFTW
using Unitful, UnitfulAstro

using AstroSimBase
using PhysicalParticles
using AstroNbodySim

using WaveDM

using CairoMakie
CairoMakie.activate!(visible = false)


# =====================================================================
# Unit tests — pure-function checks (no SPE time-stepping).
# =====================================================================

@time @testset "Unit tests" begin
    include("unit_test.jl")
end

@time @testset "MOND physics" begin
    # RAR / Milgrom / modelRAR / fit_RAR_a0 — fast (<5 s).
    include("test_MOND_units.jl")
end

@time @testset "Diagnostic helpers" begin
    # vec_cartesian_to_*, compute_madelung_velocity, compute_beta_star,
    # parse_waveDM_filename, findfirstvalue, etc.
    include("test_diagnostic_units.jl")
end

@time @testset "Catalog accessors" begin
    # Every list_*, *_index, default_RAR_radial_window branch.
    include("test_catalog_accessors.jl")
end

@time @testset "compute_timestep + compute_averaged_fields" begin
    # CFL-bound `compute_timestep` (autoset on/off), `compute_averaged_fields`
    # (average=true/false × baryon_mode=:ignored / :mesh).
    include("test_compute_timestep.jl")
end

# =====================================================================
# Integration tests — full simulate_waveDM runs at low resolution.
# =====================================================================

@time @testset "Simulation" begin
    include("integration_test.jl")
end

@time @testset "simulate_waveDM smoke (model × baryon_mode matrix)" begin
    # One smoke run per (model, baryon_mode) branch surface.  Total
    # budget <5 min on a CI box; each individual run <60 s.
    include("test_simulate_waveDM_smoke.jl")
end

@time @testset "Regression: flag_RAR_snapshot cadence" begin
    # Regression test for the `Prop.csv` last-row `RAR_a0 = NaN` bug
    # fixed in `src/simulation_main.jl` (L900-955).
    include("regression_rar_snapshot.jl")
end

# =====================================================================
# 3D SPE validation — analytic-soliton convergence tests.
# =====================================================================
#
# NOTE: The historical `run_all_3d_validation_tests(verbose=false, quick=true)`
# entry point in the original `runtests.jl` was a dead reference — the
# function is not defined in the current source tree.  The 1D/3D
# validation tests in `simulation_SPE3D_validation.jl` are 1750 lines
# of analytic-soliton checks that are run interactively, not as part
# of `Pkg.test`.  If you want them in CI, port them into proper
# `@testset` blocks first.

# @testset "3D SPE Validation Tests" begin
#     include("simulation_SPE3D_validation.jl")
#     results = run_all_3d_validation_tests(verbose=false, quick=true)
#     ...
# end

@testset "Parallel backend" begin
    # The new `ParallelBackend`/`DeviceConfig` adaptive parallel system
    # is exercised end-to-end.  The distributed Poisson test is skipped
    # when fewer than 2 workers are available.
    include("test_parallel.jl")
end
