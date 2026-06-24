# Density profile models for IC generation

"""
$(TYPEDSIGNATURES)

Sample density at a given radius using a density model.

The radius `r` is in *dimensionless code units*;
it is converted to physical lengths via `length_astro`
and the model density is converted back to code units through `density_astro`.
"""
function sampling_density(r, model, length_astro, density_astro)
    return upreferred(GalacticDynamics.density(model, r * length_astro) / density_astro)
end