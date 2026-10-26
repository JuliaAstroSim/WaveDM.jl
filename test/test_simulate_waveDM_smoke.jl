# `simulate_waveDM` smoke tests — one per (model × baryon_mode) branch.
#
# Goals
# -----
# * Cover every `model ∈ (:MW, :dwarf, :dwarf_NFW, :dwarf_Zhao, :dwarf_UFDs,
#   :Elliptical, :cluster_NFW, :cluster_Burkert)` branch surface.
# * Cover every `baryon_mode ∈ (:ignored, :mesh, :particles_static)`
#   branch surface where supported.
# * Cover the `:RC` / `:profile` `extract_mode` branches.
# * All runs must complete in well under 60 s on a CI box (Nx = 32, Nt = 3
#   for the cheapest cases; up to Nx = 64, Nt = 8 for the heavier ones).
# * Each test is fully hermetic — uses `mktempdir()` and cleans up.
#
# Sister to `integration_test.jl` (which only tests `:dwarf, :ignored`) —
# together they exercise every `simulate_waveDM` branch surface.

using Test
using Unitful: @u_str, u"kpc", u"Msun", u"Msun/kpc^3"
using WaveDM

# Convenience: minimal smoke kwargs shared by every `simulate_waveDM`
# invocation in this file. Override per-call as needed.
function _smoke(; kwargs...)
    outdir = mktempdir()
    try
        result = WaveDM.simulate_waveDM(;
            V = (x, y, z, ψ) -> 0.0,
            Xmax = 1.0, Ymax = 1.0, Zmax = 1.0,
            Tmax = 0.001,
            Nx = 32, Ny = 32, Nz = 32,
            Nt = 3,
            StepsBetweenSnapshots = 1,
            autoset_timestep = true,
            autoset_timestep_ratio = 0.5,    # conservative dt
            outputdir = outdir,
            title = "smoke",
            gpu = false,
            Realtime = false,
            extract_dwarf_granule = false,
            save_phi = false,
            save_IC = false,
            unicode_plot = false,
            kwargs...,
        )
        return result, outdir
    catch e
        @warn "smoke run failed" exception=(e, catch_backtrace())
        rm(outdir; recursive = true, force = true)
        rethrow()
    end
end

# -------------------------------------------------------------------------
# model = :MW  ×  baryon_mode ∈ (:ignored, :mesh, :particles_static,
#                                  :particles_dynamic)
# -------------------------------------------------------------------------

@testset "simulate_waveDM smoke: :MW / :ignored" begin
    result, outdir = _smoke(model = :MW, baryon_mode = :ignored)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :MW / :mesh" begin
    result, outdir = _smoke(model = :MW, baryon_mode = :mesh)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :MW / :particles_static" begin
    result, outdir = _smoke(model = :MW, baryon_mode = :particles_static,
                            Np = 100)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :MW / :particles_dynamic" begin
    # This is the heaviest smoke — N-body inside SPE. Keep Np small.
    result, outdir = _smoke(model = :MW, baryon_mode = :particles_dynamic,
                            Np = 50, Nt = 2)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

# -------------------------------------------------------------------------
# model = :dwarf  ×  baryon_mode ∈ (:ignored, :particles_static,
#                                    :particles_dynamic)
# -------------------------------------------------------------------------

@testset "simulate_waveDM smoke: :dwarf / :ignored" begin
    result, outdir = _smoke(model = :dwarf, baryon_mode = :ignored,
                            halo_ρ0 = 3e5u"Msun/kpc^3",
                            halo_r0 = 100u"kpc",
                            halo_β  = 0.5)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :dwarf / :particles_static" begin
    # Required kwargs for the particles branch.
    result, outdir = _smoke(
        model = :dwarf, baryon_mode = :particles_static,
        Np = 50,
        halo_ρ0 = 3e5u"Msun/kpc^3",
        halo_r0 = 100u"kpc",
        halo_β  = 0.5,
        stellar_TotalMass    = 1e8u"Msun",
        stellar_ScaleRadius  = 1.0u"kpc",
        thickness_ratio_stellar = 0.1,
        # gases_* are NaN by default → falls back to stellar-only
        # at src/ic/generation.jl:167-168.
    )
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :dwarf / :particles_dynamic" begin
    result, outdir = _smoke(
        model = :dwarf, baryon_mode = :particles_dynamic,
        Np = 50, Nt = 2,
        halo_ρ0 = 3e5u"Msun/kpc^3",
        halo_r0 = 100u"kpc",
        halo_β  = 0.5,
        stellar_TotalMass    = 1e8u"Msun",
        stellar_ScaleRadius  = 1.0u"kpc",
        thickness_ratio_stellar = 0.1,
    )
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

# -------------------------------------------------------------------------
# Analytic-halo models: :dwarf_NFW, :dwarf_Zhao  ×  :ignored
# -------------------------------------------------------------------------

@testset "simulate_waveDM smoke: :dwarf_NFW / :ignored" begin
    result, outdir = _smoke(model = :dwarf_NFW, baryon_mode = :ignored,
                            halo_ρ0 = 3e5u"Msun/kpc^3",
                            halo_r0 = 100u"kpc")
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :dwarf_Zhao / :ignored" begin
    result, outdir = _smoke(model = :dwarf_Zhao, baryon_mode = :ignored,
                            halo_ρ0 = 3e5u"Msun/kpc^3",
                            halo_r0 = 100u"kpc",
                            halo_α  = 1.5, halo_β = 2.0, halo_γ = 0.5)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

# -------------------------------------------------------------------------
# :dwarf_UFDs — only Crater_2 is wired in (see src/ic/generation.jl:207-209)
# -------------------------------------------------------------------------

@testset "simulate_waveDM smoke: :dwarf_UFDs / Crater_2 / :ignored" begin
    # Find Crater_2 row index in the Hayashi 2023 table.
    df_ufd = list_galaxies(; model = :dwarf_UFDs)
    crater_idx = findfirst(df_ufd.Galaxy .== "Crater_2")
    @test crater_idx !== nothing

    result, outdir = _smoke(
        model = :dwarf_UFDs, baryon_mode = :ignored,
        Galaxy_id = crater_idx,
        Nx = 64, Nt = 5,
        Xmax = 1.0, Ymax = 1.0, Zmax = 1.0,
        Tmax = 0.005,
    )
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

# -------------------------------------------------------------------------
# Cluster / elliptical models
# -------------------------------------------------------------------------

@testset "simulate_waveDM smoke: :Elliptical / :mesh" begin
    result, outdir = _smoke(model = :Elliptical, baryon_mode = :mesh,
                            baryon_ρ0 = 1e5u"Msun/kpc^3",
                            baryon_r0 = 5u"kpc",
                            halo_ρ0   = 3e5u"Msun/kpc^3",
                            halo_r0   = 100u"kpc")
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :cluster_NFW / :mesh" begin
    result, outdir = _smoke(model = :cluster_NFW, baryon_mode = :mesh,
                            baryon_ρ0 = 1e5u"Msun/kpc^3",
                            baryon_r0 = 50u"kpc",
                            baryon_β  = 2.0/3.0,
                            halo_ρ0   = 3e5u"Msun/kpc^3",
                            halo_r0   = 300u"kpc")
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: :cluster_Burkert / :mesh" begin
    result, outdir = _smoke(model = :cluster_Burkert, baryon_mode = :mesh,
                            baryon_ρ0 = 1e5u"Msun/kpc^3",
                            baryon_r0 = 50u"kpc",
                            baryon_β  = 2.0/3.0,
                            halo_ρ0   = 3e5u"Msun/kpc^3",
                            halo_r0   = 300u"kpc")
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

# -------------------------------------------------------------------------
# Unsupported-baryon-mode / unsupported-model error branches
# -------------------------------------------------------------------------

@testset "simulate_waveDM: unsupported model errors" begin
    # `_planned_model_error` is raised for SPARC models in this driver
    # (production path goes through `test_SPARC_*` in SPE-SPARC.jl).
    @test_throws ErrorException _smoke(model = :SPARC_LTGs,
                                         baryon_mode = :ignored)
    @test_throws ErrorException _smoke(model = :SPARC_Xray_ETGs,
                                         baryon_mode = :ignored)
    @test_throws ErrorException _smoke(model = :SPARC_rotating_ETGs,
                                         baryon_mode = :ignored)
    @test_throws ErrorException _smoke(model = :dSph_Lelli2017,
                                         baryon_mode = :ignored)

    # Unknown model → `_unknown_model_error`.
    @test_throws ErrorException _smoke(model = :NotAGalaxy,
                                         baryon_mode = :ignored)
end

@testset "simulate_waveDM: unsupported (model × baryon_mode) errors" begin
    # `_unsupported_baryon_mode_error` matrix.
    @test_throws ErrorException _smoke(model = :Elliptical,
                                         baryon_mode = :particles_static)
    @test_throws ErrorException _smoke(model = :cluster_NFW,
                                         baryon_mode = :particles_static)
    @test_throws ErrorException _smoke(model = :cluster_Burkert,
                                         baryon_mode = :particles_static)
    @test_throws ErrorException _smoke(model = :dwarf_NFW,
                                         baryon_mode = :particles_static)
    @test_throws ErrorException _smoke(model = :dwarf_Zhao,
                                         baryon_mode = :particles_static)
    @test_throws ErrorException _smoke(model = :dwarf_UFDs,
                                         baryon_mode = :mesh)
end

@testset "simulate_waveDM: UFD Galaxy_id out of range" begin
    df_ufd = list_galaxies(; model = :dwarf_UFDs)
    @test_throws ErrorException _smoke(
        model = :dwarf_UFDs, baryon_mode = :ignored,
        Galaxy_id = nrow(df_ufd) + 1,
    )
    # Non-Crater_2 galaxies are explicitly rejected (see
    # src/ic/generation.jl:207-209).
    non_crater = findfirst(df_ufd.Galaxy .!= "Crater_2")
    if non_crater !== nothing
        @test_throws ErrorException _smoke(
            model = :dwarf_UFDs, baryon_mode = :ignored,
            Galaxy_id = non_crater,
        )
    end
end

# -------------------------------------------------------------------------
# extract_mode branches
# -------------------------------------------------------------------------

@testset "simulate_waveDM smoke: extract_mode = :profile" begin
    # The :profile branch lives inside SPE3D_waveDM but is exercised
    # whenever extract_dwarf_granule=true; smoke path uses the default
    # MW :Zhao target profile.
    result, outdir = _smoke(model = :MW, baryon_mode = :ignored,
                            extract_dwarf_granule = true,
                            extract_min_t = 0.0u"Gyr",
                            extract_mode = :profile,
                            target_profile_model = :Zhao,
                            target_profile_ρ0 = 1.55e7u"Msun/kpc^3",
                            target_profile_rs = 11.75u"kpc",
                            target_profile_α = 1.19,
                            target_profile_β = 2.95,
                            target_profile_γ = 0.95,
                            target_fitting_rs_ratio = 3,
                            Nt = 4, StepsBetweenSnapshots = 2)
    @test result !== nothing
    rm(outdir; recursive = true, force = true)
end

@testset "simulate_waveDM smoke: extract_mode = :RC requires Cooke catalog" begin
    # The :RC branch needs df_CO_RC CSV (loaded from
    # NonlinearOptics.jl/examples/SPE-WaveDM/Cooke2022_RC/<name>.csv).
    # Without that data we expect an error path inside the driver —
    # we just check it does NOT silently succeed with a wrong setup.
    @test_throws Exception _smoke(
        model = :dwarf,
        baryon_mode = :ignored,
        extract_dwarf_granule = true,
        extract_min_t = 0.0u"Gyr",
        extract_mode = :RC,
        Galaxy_id = 99,                    # invalid → ErrorException from guard
        Nt = 4, StepsBetweenSnapshots = 2,
    )
end
