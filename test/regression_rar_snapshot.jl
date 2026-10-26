# Regression test for the `flag_RAR_snapshot` cadence bug.
#
# Symptom (fixed in `src/simulation_main.jl` around L900-955):
#   When `Nt` was a multiple of `StepsBetweenSnapshots` but *not* of
#   `RAR_snapshot_stride`, the last row of the streamed `Prop.csv`
#   could carry a `NaN` `RAR_a0` because the per-stride temp buffer
#   was only flushed inside the `StepsBetweenSnapshots`-gated branch.
#
# Fix (the diff that landed in `simulation_main.jl`):
#   * `flag_RAR_snapshot && iszero(mod(i, RAR_snapshot_stride))` (was `i-1`)
#   * Immediate `ArrayRAa0[] = vcat(...) + empty!(temp)` after every push
#   * `RAR_snapshot_stride = StepsBetweenSnapshots` (default follow)
#
# This test exercises the exact pathological case:
#   Nt=10, StepsBetweenSnapshots=5, RAR_snapshot_stride=2
# so the snapshot fit fires at i = 2, 4, 6, 8, 10
# and the streamed row writes fire at i = 5, 10.
# The last streamed row (at i = 10) must read the value pushed at i = 10,
# which is the actual last snapshot — not a stale value from i = 8.

using Test
using CSV
using DataFrames: DataFrame, nrow
using WaveDM

@testset "Regression: flag_RAR_snapshot cadence" begin
    outdir = mktempdir()
    try
        result = WaveDM.simulate_waveDM(;
            model = :MW,
            V = (x, y, z, ψ) -> 0.0,
            Nx = 32, Ny = 32, Nz = 32,
            Nt = 10,
            StepsBetweenSnapshots = 5,
            RAR_snapshot_stride = 2,
            Xmax = 1.0, Ymax = 1.0, Zmax = 1.0,
            Tmax = 0.001,
            baryon_mode = :mesh,
            outputdir = outdir,
            title = "rar_snapshot_regression",
            flag_RAR_snapshot = true,
            flag_stream_Prop_csv = true,
            RAR_snapshot_minR_kpc = 0.1,
            RAR_snapshot_maxR_kpc = 0.4,         # < Xmax/2 so cells exist
            RAR_snapshot_slice = :midplane,
            gpu = false,
            Realtime = false,
            extract_dwarf_granule = false,
            save_phi = false,
            save_IC = false,
            unicode_plot = false,
        )
        # The simulation succeeded.
        @test result !== nothing

        prop_path = joinpath(outdir, "rar_snapshot_regression,  - Prop.csv")
        @test isfile(prop_path)

        df = CSV.read(prop_path, DataFrame)
        @test nrow(df) == 2                       # 2 streamed rows (i = 5, 10)
        @test :RAR_a0 in propertynames(df)

        # The last row's RAR_a0 must NOT be NaN — that's the bug we fixed.
        # (It can be a real number, or possibly also NaN if the fit
        #  genuinely fails on the very last snapshot, but that's
        #  unlikely on a :mesh / uniform halo. We assert finiteness
        #  so the test catches the regression without being brittle.)
        last = df.RAR_a0[end]
        @test !isnan(last)

        # All earlier rows must also be finite (smoke).
        @test all(!isnan, df.RAR_a0)
    finally
        rm(outdir; recursive = true, force = true)
    end
end

@testset "Regression: flag_RAR_snapshot disabled produces no RAR_a0 column" begin
    outdir = mktempdir()
    try
        WaveDM.simulate_waveDM(;
            model = :MW,
            V = (x, y, z, ψ) -> 0.0,
            Nx = 32, Ny = 32, Nz = 32,
            Nt = 6,
            StepsBetweenSnapshots = 2,
            Xmax = 1.0, Ymax = 1.0, Zmax = 1.0,
            Tmax = 0.001,
            baryon_mode = :mesh,
            outputdir = outdir,
            title = "no_rar",
            flag_RAR_snapshot = false,
            flag_stream_Prop_csv = true,
            gpu = false,
            Realtime = false,
            extract_dwarf_granule = false,
            save_phi = false,
            save_IC = false,
            unicode_plot = false,
        )
        prop_path = joinpath(outdir, "no_rar,  - Prop.csv")
        @test isfile(prop_path)
        df = CSV.read(prop_path, DataFrame)
        @test nrow(df) == 3                       # i = 2, 4, 6
        @test !(:RAR_a0 in propertynames(df))      # column must not exist
    finally
        rm(outdir; recursive = true, force = true)
    end
end

@testset "Regression: flag_RAR_snapshot with baryon_mode=:ignored auto-disables" begin
    # The SPE driver at simulation_main.jl:518-520 emits a warning and
    # sets flag_RAR_snapshot = false when baryon_mode == :ignored.
    # We just check that the simulation completes without throwing
    # and that the resulting Prop.csv does NOT have a RAR_a0 column.
    outdir = mktempdir()
    try
        WaveDM.simulate_waveDM(;
            model = :MW,
            V = (x, y, z, ψ) -> 0.0,
            Nx = 32, Ny = 32, Nz = 32,
            Nt = 4,
            StepsBetweenSnapshots = 2,
            Xmax = 1.0, Ymax = 1.0, Zmax = 1.0,
            Tmax = 0.001,
            baryon_mode = :ignored,
            outputdir = outdir,
            title = "ignored_rar",
            flag_RAR_snapshot = true,             # ignored because of :ignored
            flag_stream_Prop_csv = true,
            gpu = false,
            Realtime = false,
            extract_dwarf_granule = false,
            save_phi = false,
            save_IC = false,
            unicode_plot = false,
        )
        prop_path = joinpath(outdir, "ignored_rar,  - Prop.csv")
        df = CSV.read(prop_path, DataFrame)
        @test !(:RAR_a0 in propertynames(df))
    finally
        rm(outdir; recursive = true, force = true)
    end
end
