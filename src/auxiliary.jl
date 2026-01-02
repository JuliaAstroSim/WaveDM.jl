"""
    optimize_massive_dwarfs_RC_stellar_fitting(r_RC, v_RC, lower, upper, ic; NumSamples_stellar, thickness_ratio_stellar)

Optimize dwarf galaxy rotation curve stellar fitting.
"""
function optimize_massive_dwarfs_RC_stellar_fitting(r_RC, v_RC, lower, upper, ic;
    NumSamples_stellar = 500,
    thickness_ratio_stellar = 0.14,
)
    function target_Stellar_RC(params;
        )
        # Sample N-body particles
        TotalMass_stellar = params[1] * 1e9u"Msun"
        ScaleRadius = params[2] * u"kpc"
        if thickness_ratio_stellar > 0
            ScaleHeight = ScaleRadius * thickness_ratio_stellar
        else
            ScaleHeight = ScaleRadius * params[3]
        end

        config_stellar = AstroIC.ExponentialDisc(;
                collection = STAR,
                NumSamples = NumSamples_stellar,
                TotalMass = TotalMass_stellar,
                ScaleRadius,
                ScaleHeight,
        )
        Random.seed!(1234)
        particles_stellar = generate(config_stellar)

        # DirectSum
        accStellar = [AstroNbodySim.compute_force_at_point(PVector(r, 0.0u"kpc", 0.0u"kpc"), particles_stellar, C.G, 1.0u"kpc") for r in r_RC]
        velStellar = ustrip.(u"km/s", sqrt.(norm.(accStellar) .* r_RC))

        constraint_RC = sum((velStellar .- ustrip.(u"km/s", v_RC)) .^ 2) / length(r_RC)
        # constraint_RC = sum(abs.(velStellar .- ustrip.(u"km/s", v_RC))) / length(r_RC)
        # @show constraint_RC
        return constraint_RC
    end

    result = Optim.optimize(
        target_Stellar_RC,
        lower, upper, ic,
        Optim.Fminbox(),
        Optim.Options(
            store_trace = true,
            # iterations = 500,
            # outer_iterations = 500,
            iterations = 200,
            outer_iterations = 200,
            # iterations = 50,
            # outer_iterations = 50,
            # x_tol = 1e-10,
            # outer_x_tol = 1e-10,
        )
    )
    trace = [state.value for state in result.trace]

    return result
end
