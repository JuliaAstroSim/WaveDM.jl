# plot_constraints_julia.jl
# ===========================
#
# Mirror of ``plot_constraints_static.py``: render the wave-CDM constraint
# catalog as a horizontal bar chart on a log-scaled x-axis.  The Julia path
# uses ``ParameterSpaceVis.plot_constraints_catalog`` directly on the
# ParameterSpace.jl DataSet (multi-interval ConstraintRow per citekey),
# demonstrating that the catalog schema is round-trippable from YAML → Julia
# → figure without going through matplotlib.
#
# Backends supported:
#   * ``CairoMakie`` — vector + raster, embeds cleanly in docs (preferred).
#   * ``GLMakie``    — interactive, useful for hover/zoom during dev.
#
# Usage
# -----
#
#     $env:JULIA_DEPOT_PATH = "<some-writable-dir>"
#     julia --project=<path-to-ParameterSpace.jl> \
#           --project=<path-to-ParameterSpaceVis.jl> \
#           docs/constraints/plot_constraints_julia.jl
#
# Inputs
# ------
# ``output/ma_constraints.pspace.yaml`` (produced by ``build_constraints_dataset.jl``).
#
# Outputs
# -------
# ``output/ma_constraints_julia.png`` / ``.pdf`` (CairoMakie backend).
# ``output/ma_constraints_julia.svg`` for docs embedding.

using ParameterSpace
using ParameterSpaceVis
using ParameterSpaceVis: plot_constraints_catalog

# Choose backend: CairoMakie for static output, GLMakie for interactive.
using CairoMakie

const HERE     = abspath(@__DIR__)
const IN_YAML  = joinpath(HERE, "output", "ma_constraints.pspace.yaml")
const OUT_PNG  = joinpath(HERE, "output", "ma_constraints_julia.png")
const OUT_PDF  = joinpath(HERE, "output", "ma_constraints_julia.pdf")
const OUT_SVG  = joinpath(HERE, "output", "ma_constraints_julia.svg")

function main()
    isfile(IN_YAML) || error("Missing $IN_YAML; run build_constraints_dataset.jl first.")
    @info "Loading DataSet" path=IN_YAML
    ds = load_dataset(IN_YAML)

    # Filter to x-axis range (drop entries whose intervals don't intersect).
    # In this catalog every row is in-range so we keep everything; if needed
    # we could `filter(r -> any(iv -> iv.upper >= -26 && iv.lower <= -18, r.intervals), ds.data)`.

    @info "Plotting catalog" rows=length(ds)
    fig = plot_constraints_catalog(ds;
        param     = :log10_m_a_eV,
        x_min     = -26.0,
        x_max     = -18.0,
        row_label_fontsize = 10.0,
        bar_height         = 0.40,
        legend_position    = :top,
    )

    @info "Saving" OUT_PNG OUT_PDF OUT_SVG
    save(OUT_PNG, fig; pt_per_unit = 2)
    save(OUT_PDF, fig)
    save(OUT_SVG, fig)

    @info "Done."
    return fig
end

main()