# Catalog / index / default-window unit tests.
#
# Every public catalog accessor (`list_*`, `*_index`,
# `default_RAR_radial_window`) gets a smoke + out-of-range test here so
# that any breakage in `src/ic/catalog.jl` is caught immediately.
#
# Sister to `unit_test.jl` (which only covers the Lelli 2017 catalog) —
# together they exercise every catalog branch surface.

using Test
using WaveDM: list_supported_models, list_galaxies, print_catalog,
               list_MW_satellites, list_SPARC_LTGs,
               list_SPARC_Xray_ETGs, list_SPARC_rotating_ETGs,
               list_cooke_dwarfs, list_Lelli2017_dsphs,
               ufd_index, MW_satellite_index,
               SPARC_LTG_index, SPARC_Xray_ETG_index,
               SPARC_rotating_ETG_index, cooke_dwarf_index,
               Lelli2017_dsph_index,
               default_RAR_radial_window, _resolve_rar_window
using DataFrames: DataFrame, nrow

@testset "Catalog: list_supported_models" begin
    models = list_supported_models()
    @test models isa Vector
    @test !isempty(models)
    syms = [m.model for m in models]
    # All currently-implemented branches must be present
    for m in (:MW, :dwarf, :dwarf_NFW, :dwarf_Zhao, :dwarf_UFDs,
              :Elliptical, :cluster_NFW, :cluster_Burkert)
        @test m in syms
    end
    # `:MW` should advertise at least one allowed `baryon_mode`
    mw = findfirst(m -> m.model == :MW, models)
    @test mw !== nothing
    @test :ignored in models[mw].baryon_mode
    @test :particles_static in models[mw].baryon_mode
end

@testset "Catalog: list_galaxies(model = :MW)" begin
    df = list_galaxies(; model = :MW)
    @test df isa DataFrame
    @test nrow(df) >= 1
end

@testset "Catalog: list_MW_satellites + MW_satellite_index" begin
    df = list_MW_satellites()
    @test df isa DataFrame
    @test nrow(df) >= 1

    # In-range id resolves to the matching row
    @test MW_satellite_index(df, 1) == 1
    @test MW_satellite_index(df, nrow(df)) == nrow(df)

    # Out-of-range ids throw the standard error
    @test_throws ErrorException MW_satellite_index(df, 0)
    @test_throws ErrorException MW_satellite_index(df, nrow(df) + 1)
    @test_throws ErrorException MW_satellite_index(df, -1)
end

@testset "Catalog: list_SPARC_LTGs + SPARC_LTG_index" begin
    df = list_SPARC_LTGs()
    @test df isa DataFrame
    @test nrow(df) >= 1
    @test SPARC_LTG_index(df, 1) == 1
    @test SPARC_LTG_index(df, nrow(df)) == nrow(df)
    @test_throws ErrorException SPARC_LTG_index(df, 0)
    @test_throws ErrorException SPARC_LTG_index(df, nrow(df) + 1)
end

@testset "Catalog: list_SPARC_Xray_ETGs + SPARC_Xray_ETG_index" begin
    df = list_SPARC_Xray_ETGs()
    @test df isa DataFrame
    @test nrow(df) >= 1
    @test SPARC_Xray_ETG_index(df, 1) == 1
    @test SPARC_Xray_ETG_index(df, nrow(df)) == nrow(df)
    @test_throws ErrorException SPARC_Xray_ETG_index(df, 0)
    @test_throws ErrorException SPARC_Xray_ETG_index(df, nrow(df) + 1)
end

@testset "Catalog: list_SPARC_rotating_ETGs + SPARC_rotating_ETG_index" begin
    df = list_SPARC_rotating_ETGs()
    @test df isa DataFrame
    @test nrow(df) >= 1
    @test SPARC_rotating_ETG_index(df, 1) == 1
    @test SPARC_rotating_ETG_index(df, nrow(df)) == nrow(df)
    @test_throws ErrorException SPARC_rotating_ETG_index(df, 0)
    @test_throws ErrorException SPARC_rotating_ETG_index(df, nrow(df) + 1)
end

@testset "Catalog: list_cooke_dwarfs + cooke_dwarf_index" begin
    df = list_cooke_dwarfs()
    @test df isa DataFrame
    @test nrow(df) >= 1
    # Cooke dwarfs are looked up by *name* (string), not integer id.
    # Take the first row's Galaxy and round-trip it.
    name = df.Galaxy[1]
    @test cooke_dwarf_index(name) == 1
    @test_throws ErrorException cooke_dwarf_index("NotARealGalaxyName")
end

@testset "Catalog: ufd_index (Hayashi 2023 UFDs)" begin
    df = list_galaxies(; model = :dwarf_UFDs)
    @test df isa DataFrame
    @test nrow(df) >= 1
    @test ufd_index(df, 1) == 1
    @test ufd_index(df, nrow(df)) == nrow(df)
    @test_throws ErrorException ufd_index(df, 0)
    @test_throws ErrorException ufd_index(df, nrow(df) + 1)
end

@testset "default_RAR_radial_window: catalog branches" begin
    # :MW — falls back to the user-supplied floor/ceiling (no catalog default).
    rmin, rmax = default_RAR_radial_window(:MW)
    @test rmin == 0.0
    @test rmax == 0.0
    # Explicit floor / ceiling is honored.
    rmin2, rmax2 = default_RAR_radial_window(:MW; minR_kpc = 5.0, maxR_kpc = 25.0)
    @test rmin2 == 5.0
    @test rmax2 == 25.0

    # :dwarf_UFDs — Hayashi 2023 `b`: rMin = b/2, rMax = 4b.
    df_ufd = list_galaxies(; model = :dwarf_UFDs)
    i = 1
    b_kpc = ustrip(u"kpc", df_ufd.b[i])
    rmin, rmax = default_RAR_radial_window(:dwarf_UFDs; Galaxy_id = i)
    @test rmin ≈ b_kpc / 2
    @test rmax ≈ 4 * b_kpc
    # Out-of-range id → ArgumentError
    @test_throws ArgumentError default_RAR_radial_window(:dwarf_UFDs;
                                                          Galaxy_id = nrow(df_ufd) + 1)

    # :SPARC_LTGs — rMax = 2 × max(rs, RHI)
    df_ltg = list_SPARC_LTGs()
    rmin, rmax = default_RAR_radial_window(:SPARC_LTGs; Galaxy_id = 1)
    rs_kpc  = Float64(df_ltg.rs[1])
    RHI_kpc = Float64(df_ltg.RHI[1])
    @test rmin == 0.0
    @test rmax ≈ 2 * max(rs_kpc, RHI_kpc)
    @test_throws ArgumentError default_RAR_radial_window(:SPARC_LTGs;
                                                          Galaxy_id = nrow(df_ltg) + 1)

    # :SPARC_Xray_ETGs — rMin = 0.1 Reff, rMax = 4 Reff
    df_xray = list_SPARC_Xray_ETGs()
    rmin, rmax = default_RAR_radial_window(:SPARC_Xray_ETGs; Galaxy_id = 1)
    Reff_kpc = Float64(df_xray.Reff[1])
    @test rmin ≈ 0.1 * Reff_kpc
    @test rmax ≈ 4.0 * Reff_kpc

    # :SPARC_rotating_ETGs — rMax = 4 × max(Reff, Rexp)
    df_rot = list_SPARC_rotating_ETGs()
    rmin, rmax = default_RAR_radial_window(:SPARC_rotating_ETGs; Galaxy_id = 1)
    Reff_kpc = Float64(df_rot.Reff[1])
    Rexp_kpc = Float64(df_rot.Rexp[1])
    @test rmin == 0.0
    @test rmax ≈ 4 * max(Reff_kpc, Rexp_kpc)

    # Unknown catalog — returns the user-supplied floor/ceiling unchanged.
    rmin, rmax = default_RAR_radial_window(:NotARealCatalog;
                                             minR_kpc = 3.0, maxR_kpc = 12.0)
    @test rmin == 3.0
    @test rmax == 12.0
end

@testset "_resolve_rar_window: :auto sentinel" begin
    # With :auto the per-catalog default kicks in for that bound; the
    # other bound (numeric) acts as the floor / ceiling.
    df_ltg = list_SPARC_LTGs()
    rs_kpc  = Float64(df_ltg.rs[1])
    RHI_kpc = Float64(df_ltg.RHI[1])
    expected_rmax = 2 * max(rs_kpc, RHI_kpc)

    # minR_kpc = :auto, maxR_kpc = 30 (numeric ceiling)
    rmin, rmax = _resolve_rar_window(:auto, 30.0, :SPARC_LTGs; Galaxy_id = 1)
    @test rmin == 0.0                                                # SPARC default
    @test rmax ≈ min(expected_rmax, 30.0)                            # ceiling applies

    # minR_kpc = 5 (numeric floor), maxR_kpc = :auto
    rmin2, rmax2 = _resolve_rar_window(5.0, :auto, :SPARC_LTGs; Galaxy_id = 1)
    @test rmin2 == 5.0                                               # floor raises rMin
    @test rmax2 ≈ expected_rmax                                      # catalog default

    # Both :auto → catalog default, no floor/ceiling.
    rmin3, rmax3 = _resolve_rar_window(:auto, :auto, :SPARC_LTGs; Galaxy_id = 1)
    @test rmin3 == 0.0
    @test rmax3 ≈ expected_rmax
end

@testset "print_catalog: smoke" begin
    # Just make sure it doesn't throw; the function prints to stdout.
    @test nothing === print_catalog()
end
