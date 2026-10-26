# Diagnostic / utility unit tests.
#
# Covers the smaller exports that don't need a full SPE run:
# * `parse_waveDM_filename`, `extract_ma_number`   (src/utils/)
# * `vec_cartesian_to_spherical`, `vec_cartesian_to_cylindrical`
# * `compute_madelung_velocity`, `compute_kinetic_bin_profiles`
# * `compute_beta_star`, `compute_mass_weighted_center`,
#   `compute_interior_mass`
# * `findfirstvalue`
#
# These are the testset gaps called out in the test-coverage survey —
# none of them is exercised by `unit_test.jl`.

using Test
using Statistics
using WaveDM: parse_waveDM_filename, extract_ma_number,
               vec_cartesian_to_spherical, vec_cartesian_to_cylindrical,
               compute_madelung_velocity, compute_kinetic_bin_profiles,
               compute_beta_star, compute_mass_weighted_center,
               compute_interior_mass, findfirstvalue

# -------------------------------------------------------------------------
# Filename parsing
# -------------------------------------------------------------------------

@testset "parse_waveDM_filename: full" begin
    fname = "milky_way_ma(1.50)_Nx(256)_Nt(1000)_Xmax(20.0).jld2"
    parsed = parse_waveDM_filename(fname)
    @test parsed.Nx == 256
    @test parsed.Nt == 1000
    @test parsed.ma == 1.5
    @test parsed.Xmax == 20.0
end

@testset "parse_waveDM_filename: missing tags" begin
    # Any absent tag becomes nothing.
    parsed = parse_waveDM_filename("foo_ma(2.0)_Nx(128).jld2")
    @test parsed.Nx == 128
    @test parsed.ma == 2.0
    @test parsed.Nt === nothing
    @test parsed.Xmax === nothing

    parsed2 = parse_waveDM_filename("no_tags_here.jld2")
    @test parsed2.Nx === nothing
    @test parsed2.Nt === nothing
    @test parsed2.ma === nothing
    @test parsed2.Xmax === nothing
end

@testset "extract_ma_number" begin
    @test extract_ma_number("foo_ma(3.14)_bar.jld2") == 3.14
    @test extract_ma_number("ma(0.5)") == 0.5
    @test extract_ma_number("no_ma_tag") === nothing
    @test extract_ma_number("") === nothing
end

# -------------------------------------------------------------------------
# Vector field conversions
# -------------------------------------------------------------------------

@testset "vec_cartesian_to_spherical: pure-radial field" begin
    # On a pure radial field (ax = r, ay = az = 0), the radial component
    # should match sqrt(ax^2 + ay^2 + az^2) at every cell.
    N = 8
    x = collect(LinRange(-1.0, 1.0, N))
    y = collect(LinRange(-1.0, 1.0, N))
    z = collect(LinRange(-1.0, 1.0, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)
    rrr = sqrt.(xxx.^2 .+ yyy.^2 .+ zzz.^2)

    ax = xxx; ay = zeros(N, N, N); az = zeros(N, N, N)

    ar, aθ, aϕ = vec_cartesian_to_spherical(ax, ay, az, xxx, yyy, zzz, rrr)
    # ar should equal |v| at every cell.
    @test all(ar .≈ sqrt.(ax.^2 .+ ay.^2 .+ az.^2))

    # Polar/azimuthal components vanish for a pure radial field.
    @test all(abs.(aθ) .< 1e-10)
    @test all(abs.(aϕ) .< 1e-10)
end

@testset "vec_cartesian_to_cylindrical: pure-azimuthal field" begin
    # v = (-y, x, 0) is the canonical azimuthal unit vector (up to /r).
    N = 8
    x = collect(LinRange(-1.0, 1.0, N))
    y = collect(LinRange(-1.0, 1.0, N))
    z = collect(LinRange(-1.0, 1.0, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)
    rrr_xy = sqrt.(xxx.^2 .+ yyy.^2)

    ax = -yyy; ay = xxx; az = zeros(N, N, N)

    aρ, aϕ, az_cyl = vec_cartesian_to_cylindrical(ax, ay, az, xxx, yyy)
    # The radial component must vanish (v ⊥ ρ̂).
    @test all(abs.(aρ) .< 1e-10)
    # |v|² = x² + y² = rxy², and aϕ = |v| at every cell.
    @test all(abs.(aϕ) .- sqrt.(ax.^2 .+ ay.^2) .< 1e-10)
    # z component passes through unchanged.
    @test all(az_cyl .≈ az)
end

# -------------------------------------------------------------------------
# Madelung velocity reconstruction
# -------------------------------------------------------------------------

@testset "compute_madelung_velocity: known field" begin
    # Pick a wavefunction whose phase gradient is uniform: ψ(r) =
    # exp(i k·r) * Gaussian envelope.  Then v = ∇phase = k (in the bulk).
    N = 16
    L = 4.0
    x = collect(LinRange(-L/2, L/2, N))
    y = collect(LinRange(-L/2, L/2, N))
    z = collect(LinRange(-L/2, L/2, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)

    kvec = (0.7, -0.4, 0.2)            # code-unit velocity = k
    env  = exp.(-(xxx.^2 .+ yyy.^2 .+ zzz.^2) ./ 2)
    ψ = env .* exp.(im .* (kvec[1] .* xxx .+ kvec[2] .* yyy .+ kvec[3] .* zzz))

    vx, vy, vz = compute_madelung_velocity(ψ, xxx, yyy, zzz)

    # Check the central bulk (away from the Gaussian edges): each v_? ≈ k_?.
    # We use the central 4⁰-cell cube (8 kpc³ region).
    N2 = 6
    mid = div(N, 2)
    sl = (mid - N2):(mid + N2 - 1)
    @test all(isapprox.(vx[sl, sl, sl], kvec[1]; atol = 0.3))
    @test all(isapprox.(vy[sl, sl, sl], kvec[2]; atol = 0.3))
    @test all(isapprox.(vz[sl, sl, sl], kvec[3]; atol = 0.3))

    # All velocities must be finite.
    @test all(isfinite, vx); @test all(isfinite, vy); @test all(isfinite, vz)
end

# -------------------------------------------------------------------------
# Mass center / interior mass
# -------------------------------------------------------------------------

@testset "compute_mass_weighted_center: shifted Gaussian" begin
    # A Gaussian density centered at (1, -2, 0.5) should report those
    # coordinates (in code units) as its mass-weighted center.
    N = 16
    L = 4.0
    x = collect(LinRange(-L/2, L/2, N))
    y = collect(LinRange(-L/2, L/2, N))
    z = collect(LinRange(-L/2, L/2, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)
    Δx = x[2] - x[1]
    r2 = (xxx .- 1.0).^2 .+ (yyy .+ 2.0).^2 .+ (zzz .- 0.5).^2
    ρ = exp.(-r2 ./ 0.5)

    cx, cy, cz = compute_mass_weighted_center(ρ, xxx, yyy, zzz, Δx)
    @test cx ≈ 1.0  atol = 1e-2
    @test cy ≈ -2.0 atol = 1e-2
    @test cz ≈ 0.5  atol = 1e-2
end

@testset "compute_interior_mass: enclosed fraction of a Gaussian" begin
    N = 32
    L = 8.0
    x = collect(LinRange(-L/2, L/2, N))
    y = collect(LinRange(-L/2, L/2, N))
    z = collect(LinRange(-L/2, L/2, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)
    rrr = sqrt.(xxx.^2 .+ yyy.^2 .+ zzz.^2)
    Δx = x[2] - x[1]
    σ = 0.5
    ρ = exp.(-rrr.^2 ./ (2σ^2))

    # Mass inside r < σ should be ~ 1 - exp(-1/2) ≈ 0.393 of the total.
    M_total = sum(ρ) * Δx^3
    M_inside = compute_interior_mass(ρ, rrr, Δx^3, 1.0; r_max = σ)
    @test M_inside / M_total ≈ 1 - exp(-1/2)  rtol = 5e-3
end

# -------------------------------------------------------------------------
# Inner-slope fit (compute_beta_star)
# -------------------------------------------------------------------------

@testset "compute_beta_star: power-law density" begin
    # ρ(r) ∝ r^α → log ρ = α log r + c → fit slope = α.
    N = 32
    L = 8.0
    x = collect(LinRange(-L/2, L/2, N))
    y = collect(LinRange(-L/2, L/2, N))
    z = collect(LinRange(-L/2, L/2, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)
    rrr = sqrt.(xxx.^2 .+ yyy.^2 .+ zzz.^2)

    for α_test in (-1.5, -1.0, -0.5)
        ρ = rrr.^α_test
        β_fit = compute_beta_star(rrr, ρ, 0.5u"kpc", 2.0u"kpc", 1.0u"kpc")
        @test β_fit ≈ α_test  rtol = 5e-2
    end
end

# -------------------------------------------------------------------------
# Kinetic-energy-bin profiles
# -------------------------------------------------------------------------

@testset "compute_kinetic_bin_profiles: single-bin mass fraction" begin
    # All mass at v = 5 km/s inside r_max should land in the v_edges bin
    # that contains 5 km/s.  With v_edges = [0, 10], that bin is the
    # first one → mass_fraction_inner + mass_fraction_annulus = 1.
    N = 8
    L = 2.0
    x = collect(LinRange(-L/2, L/2, N))
    y = collect(LinRange(-L/2, L/2, N))
    z = collect(LinRange(-L/2, L/2, N))
    oneM = ones(N, N, N)
    xxx = oneM .* x; yyy = oneM .* y'; zzz = oneM .* reshape(z, 1, 1, N)
    Δx = x[2] - x[1]
    rrr = sqrt.(xxx.^2 .+ yyy.^2 .+ zzz.^2)
    ρ = ones(N, N, N) .* 1.0               # uniform density
    vx = ones(N, N, N) .* 5.0              # all moving at 5 km/s in x
    vy = zeros(N, N, N)
    vz = zeros(N, N, N)

    result = compute_kinetic_bin_profiles(rrr, ρ, vx, vy, vz;
        r_max = 0.5u"kpc", v_bin_edges = [0.0, 10.0, 20.0],
        section = 1, uL = 1.0, uV = 1.0, uDensity = 1.0,
        r_min_mass = 0.1u"kpc", r_max_mass = 0.4u"kpc",
    )

    @test haskey(result, :mass_fraction_inner)
    @test haskey(result, :mass_fraction_annulus)
    @test result.mass_fraction_inner[1] + result.mass_fraction_annulus[1] ≈ 1.0
    # The second bin (10–20 km/s) should have zero mass.
    @test result.mass_fraction_annulus[2] == 0.0
end

# -------------------------------------------------------------------------
# findfirstvalue (src/simulation/utils.jl re-export)
# -------------------------------------------------------------------------

@testset "findfirstvalue: matching index" begin
    @test findfirstvalue([10.0, 20.0, 30.0], 20.0) == 2
    @test findfirstvalue([10.0, 20.0, 30.0], 30.0) == 3
    @test findfirstvalue([10.0, 20.0, 30.0], 10.0) == 1
end

@testset "findfirstvalue: no match returns nothing" begin
    @test findfirstvalue([10.0, 20.0, 30.0], 25.0) === nothing
    @test findfirstvalue(Float64[], 5.0) === nothing
end
