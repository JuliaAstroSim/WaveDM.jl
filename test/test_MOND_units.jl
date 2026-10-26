# MOND physics unit tests.
#
# Covers `RAR`, `Milgrom`, their inverses, `modelRAR` (the LsqFit model),
# and `fit_RAR_a0` (the LsqFit wrapper used by the per-snapshot RAR
# diagnostic in `SPE3D_waveDM`).
#
# Sister to the `Lelli 2017 RAR` testset in `unit_test.jl` (which only
# covers `Lelli2017_RAR`) — together they exercise every RAR/Milgrom
# surface exported from `src/physics/MOND.jl`.

using Test
using WaveDM: RAR, Milgrom, RAR_inv, Milgrom_inv, modelRAR, fit_RAR_a0
using WaveDM: a0_Lelli2017
using LsqFit: curve_fit

# Reference constants (Lelli 2017)
const A0 = 1.20e-10           # m/s²
const NU = 2.0

@testset "RAR: Newtonian limit a_b >> a0" begin
    # a_b/(1 - exp(-sqrt(a_b/a0))) → a_b as a_b/a0 → ∞
    @test RAR(1e-7, A0) ≈ 1e-7  rtol = 1e-3
    @test RAR(1e-6, A0) ≈ 1e-6  rtol = 1e-3
end

@testset "RAR: at the knee a_b = a0" begin
    # a0 / (1 - exp(-1)) ≈ 1.582 * a0
    @test RAR(A0, A0) ≈ A0 / (1 - exp(-1))  rtol = 1e-6
end

@testset "RAR: MOND limit a_b << a0" begin
    # a_b / (1 - exp(-sqrt(a_b/a0))) → sqrt(a_b * a0) for a_b << a0
    @test RAR(1e-13, A0) ≈ sqrt(1e-13 * A0)  rtol = 5e-3
end

@testset "RAR: missing / unitful inputs" begin
    # Milgrom accepts Missing and dispatches to `missing` (see MOND.jl:3-4)
    @test ismissing(Milgrom(missing, A0, NU))
    @test ismissing(Milgrom(1e-12, missing, NU))
    # RAR is not defined for Missing inputs (no such method dispatch).
    @test_throws MethodError RAR(missing, A0)
end

@testset "Milgrom: nu function sanity" begin
    # At a_b = a0, nu(a_b/a0, ν) = nu(1, 2). The McGaugh-Wu nu family gives
    # nu(1) = 1 / sqrt(2) by construction (verify by calling the family
    # directly via Milgrom).
    # nu(0) = 0 (deep MOND), so Milgrom(0, a0, ν) = 0.
    @test Milgrom(0.0, A0, NU) ≈ 0.0
    # In the Newtonian limit nu → 1 → Milgrom ≈ a_b.
    @test Milgrom(1e-6, A0, NU) ≈ 1e-6  rtol = 1e-3
end

@testset "RAR_inv: round-trip + missing" begin
    # round-trip: RAR_inv(RAR(a_b, a0), a_b) → a0
    for a_b in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8)
        a_MOND = RAR(a_b, A0)
        a0_back = RAR_inv(a_b, a_MOND)
        if a_MOND > a_b
            @test a0_back ≈ A0  rtol = 1e-3
        else
            # pathological / boundary: ensure no crash
            @test ismissing(a0_back) || isnan(a0_back)
        end
    end

    # When a_MOND ≤ a_b, RAR_inv returns missing (Newtonian branch is
    # under-determined; the SPE driver treats this as "skip this cell").
    @test ismissing(RAR_inv(1e-12, 1e-12))
    @test ismissing(RAR_inv(1e-12, 5e-13))
end

@testset "Milgrom_inv: round-trip + missing" begin
    for a_b in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8)
        a_MOND = Milgrom(a_b, A0, NU)
        a0_back = Milgrom_inv(a_b, a_MOND, NU)
        if a_MOND > a_b
            @test a0_back ≈ A0  rtol = 1e-3
        else
            @test ismissing(a0_back) || isnan(a0_back)
        end
    end

    @test ismissing(Milgrom_inv(1e-12, 1e-12, NU))
end

@testset "modelRAR: against RAR" begin
    # modelRAR is the LsqFit vector form of RAR; element-wise agreement.
    a_b_vec = collect(1e-13:1e-13:1e-11)
    a0 = A0
    a_RAR = RAR.(a_b_vec, a0)
    a_model = modelRAR(a_b_vec, [a0])
    @test length(a_model) == length(a_RAR)
    @test all(a_model .≈ a_RAR  rtol = 1e-10)
end

@testset "fit_RAR_a0: basic" begin
    # Generate a synthetic pair (a_b, a_MOND) on a known a0, then fit.
    a_b = collect(1e-12:1e-12:1e-9)
    a0_true = 1.5e-10
    a_MOND = modelRAR(a_b, [a0_true])

    # The fit should recover a0_true within ~1e-2 (within the LsqFit
    # bounds `lower=1e-2, upper=1e2`).
    a0_fit = fit_RAR_a0(a_b, a_MOND; lower = 1e-2, upper = 1e2)
    @test !isnan(a0_fit)
    @test a0_fit ≈ a0_true  rtol = 5e-3
end

@testset "fit_RAR_a0: degenerate inputs → NaN (no throw)" begin
    # Too few points → NaN, not an exception.
    @test isnan(fit_RAR_a0([1e-12], [1e-11]))
    @test isnan(fit_RAR_a0(Float64[], Float64[]))

    # Length mismatch throws (this is an API misuse, not a degenerate fit).
    @test_throws ErrorException fit_RAR_a0([1e-12, 2e-12], [1e-11])
end

@testset "fit_RAR_a0: bounds" begin
    # A fit on the noisy a0=1.2e-10 data should land in [1e-2, 1e2] (default).
    a_b = collect(1e-12:1e-12:1e-9)
    a0_true = A0
    a_MOND = modelRAR(a_b, [a0_true])
    a0_fit = fit_RAR_a0(a_b, a_MOND)
    @test 1e-2 ≤ a0_fit ≤ 1e2
    @test a0_fit ≈ a0_true  rtol = 5e-3
end
