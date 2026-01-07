@testset "MWE: Dwarf galaxy" begin
    astro()
    C = Constant(uAstro)
    
    model = :dwarf
    V = (x,y,z,ψ)->0.0
    FDM_mass_ratio = 1.0
    FDM_radius_ratio = 1.0
    title = "Test_Dwarf_Minimal"
    Xmax = 1.0
    Ymax = 1.0
    Zmax = 1.0
    Tmax = 0.001
    Nt = 3
    Nx = 32
    Ny = Nx
    Nz = Nx
    Np = 100
    absorb_coeff = 10
    StepsBetweenSnapshots = 1
    IC_vel = nothing
    reset_velocity = false
    Φ_b = nothing
    ax_b = nothing
    ay_b = nothing
    az_b = nothing
    IC_only = false
    static = false
    save_IC = false
    rotational_ratio = 0.0
    velocity_ratio = 1.0
    velocity_falling = false
    outputdir = joinpath(tempdir(), "WaveDM_test")
    massRadius = 10.0u"kpc"
    bulk_perturb = true
    bulk_size = 4
    bulk_center_size = 0
    bulk_shift_size = 1
    baryon_β = 0.5
    baryon_ρ0 = 3e5u"Msun/kpc^3"
    baryon_r0 = 100u"kpc"
    halo_β = 0.5
    halo_ρ0 = 3e5u"Msun/kpc^3"
    halo_r0 = 100u"kpc"
    halo_α = 1.0
    halo_γ = 1.0
    halo_Q = 1.0
    baryon_fraction_limit = 1
    GravitySolver = Tree()
    boundary = Periodic()
    SofteningLength = 1.0u"kpc"
    baryon_mode = :ignored
    gpu = false
    
    stellar_TotalMass = NaN
    stellar_ScaleRadius = NaN
    thickness_ratio_stellar = NaN
    gases_TotalMass = NaN
    gases_ScaleRadius = NaN
    thickness_ratio_gases = NaN
    
    mₐ = 0.2 * 1e-22 * 1.783e-36u"kg"
    Ωₘ₀ = 0.31
    aₛ = +8.29e-60u"fm"
    
    result = test_MW_MOND(;
        model, V, FDM_mass_ratio, FDM_radius_ratio, title,
        Xmax, Ymax, Zmax, Tmax, Nt, Nx, Ny, Nz, Np,
        absorb_coeff, StepsBetweenSnapshots, IC_vel, reset_velocity,
        Φ_b, ax_b, ay_b, az_b, IC_only, static, save_IC,
        rotational_ratio, velocity_ratio, velocity_falling,
        outputdir, massRadius, bulk_perturb, bulk_size,
        bulk_center_size, bulk_shift_size,
        baryon_β, baryon_ρ0, baryon_r0,
        halo_β, halo_ρ0, halo_r0, halo_α, halo_γ, halo_Q,
        baryon_fraction_limit, GravitySolver, boundary, SofteningLength,
        baryon_mode, gpu,
        stellar_TotalMass, stellar_ScaleRadius, thickness_ratio_stellar,
        gases_TotalMass, gases_ScaleRadius, thickness_ratio_gases,
        mₐ, Ωₘ₀, aₛ,
        distributed_memory = false,
        Realtime = true,
        unicode_plot = false,
        extract_dwarf_granule = false,
        save_phi = true,
        extract_min_t = 1.5u"Gyr",
        target_profile_model = :dwarf_gNFW,
        target_profile_β = halo_β,
        target_profile_β_u = halo_β + 0.1,
        target_profile_β_d = halo_β - 0.1,
        target_profile_ρ0 = halo_ρ0,
        target_profile_ρ0_u = halo_ρ0 * 1.1,
        target_profile_ρ0_d = halo_ρ0 * 0.9,
        target_profile_rs = halo_r0,
        target_profile_rs_u = halo_r0 * 1.1,
        target_profile_rs_d = halo_r0 * 0.9,
        target_fitting_rs_ratio = 3,
        target_profile_error = false,
        extract_mode = :profile,
        target_beta_star = NaN,
        target_beta_star_u = NaN,
        target_beta_star_d = NaN,
        beta_star_error_threshold = NaN,
        plot_virial = false,
    )
    
    @test !isnothing(result)
end
