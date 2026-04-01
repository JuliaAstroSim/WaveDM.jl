# Density profile models for IC generation

"""
$(TYPEDSIGNATURES)

Sample density at a given radius using a density model.
"""
function sampling_density(r, model, length_astro, density_astro)
    out_density = upreferred(GalacticDynamics.density(model, r*length_astro) / density_astro)
end