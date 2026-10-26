#!/usr/bin/env julia
# Direct runner for the new WaveDM.jl CI tests.
#
# Usage (from `E:\JuliaAstroSim`):
#   $env:JULIA_DEPOT_PATH = "<your-depot>"
#   julia --project=WaveDM.jl WaveDM.jl/test/run_all_ci.jl
#
# Or simply use `Pkg.test("WaveDM")` — both will run the same testsets.
#
# This script runs each testset in isolation with timing so that a
# partial failure (e.g. one broken module) doesn't block the rest.

using Test
using WaveDM
using CairoMakie
CairoMakie.activate!(visible = false)

# Pretty banner between testsets so output is grep-friendly.
function banner(name)
    println()
    println("=" ^ 70)
    println("  TEST SET: ", name)
    println("=" ^ 70)
end

# Track cumulative time.
const t_total = @elapsed begin
end
total_sec = Ref(0.0)

function run(name::AbstractString; include_path::AbstractString)
    banner(name)
    t0 = time()
    @testset "$name" begin
        include(include_path)
    end
    dt = time() - t0
    total_sec[] += dt
    println("✓ $name finished in $(round(dt; digits=1))s")
end

# =====================================================================
# UNIT-LEVEL (no full SPE time-stepping)
# =====================================================================

run("Unit tests (existing)"; include_path = joinpath(@__DIR__, "unit_test.jl"))
run("MOND physics";          include_path = joinpath(@__DIR__, "test_MOND_units.jl"))
run("Diagnostic helpers";     include_path = joinpath(@__DIR__, "test_diagnostic_units.jl"))
run("Catalog accessors";      include_path = joinpath(@__DIR__, "test_catalog_accessors.jl"))
run("compute_timestep";       include_path = joinpath(@__DIR__, "test_compute_timestep.jl"))

# =====================================================================
# INTEGRATION (full simulate_waveDM runs at low resolution)
# =====================================================================

run("Integration: dwarf :ignored (existing)"; include_path = joinpath(@__DIR__, "integration_test.jl"))
run("Smoke: model × baryon_mode matrix";     include_path = joinpath(@__DIR__, "test_simulate_waveDM_smoke.jl"))
run("Regression: flag_RAR_snapshot cadence"; include_path = joinpath(@__DIR__, "regression_rar_snapshot.jl"))

# =====================================================================
# PARALLEL BACKEND
# =====================================================================

run("Parallel backend"; include_path = joinpath(@__DIR__, "test_parallel.jl"))

banner("SUMMARY")
println("Total wallclock: $(round(total_sec[]; digits=1))s")
println("All done.")
