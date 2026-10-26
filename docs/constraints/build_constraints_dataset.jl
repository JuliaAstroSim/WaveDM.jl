# build_constraints_dataset.jl
# ============================
#
# Reads ``ma_constraints.yaml`` (produced by ``ma_constraints_data.py``) and
# rebuilds the constraint set using **ParameterSpace.jl**'s ``ConstraintRow``
# + ``DataSet`` API (the catalog-friendly multi-interval type).  The output is
# a second YAML file (``ma_constraints.pspace.yaml``) that follows
# ParameterSpace's serialisation format.
#
# Each YAML row (= one citekey) becomes ONE ``ConstraintRow`` carrying ALL of
# its intervals (a single paper can constrain both ``m < a`` and ``m > b``),
# plus the free-text ``data`` and ``method`` strings.  Authors / DOI / eprint
# / bibcode / url from the BibTeX parsing are propagated onto the
# ``DataSource``.
#
# Usage
# -----
#
#     $env:JULIA_DEPOT_PATH = "<some-writable-dir>"
#     julia --project=<path-to-ParameterSpace.jl> docs/constraints/build_constraints_dataset.jl
#
# The script writes:
#
# * ``output/ma_constraints.pspace.yaml`` — ParameterSpace-style YAML using
#   ``ConstraintRow`` (catalog schema).
# * ``output/ma_constraints_filtered.yaml`` — example of filtering the catalog
#   by year range.

using ParameterSpace
using ParameterSpace: ConstraintInterval, ConstraintRow, DataSource, DataSet,
                      is_include, is_exclude, save_dataset, add!
using YAML

const HERE     = abspath(@__DIR__)
const IN_YAML  = joinpath(HERE, "output", "ma_constraints.yaml")
const OUT_YAML = joinpath(HERE, "output", "ma_constraints.pspace.yaml")
const FILTERED = joinpath(HERE, "output", "ma_constraints_filtered.yaml")

"""
    _interval_kind(yaml_type::String) :: Symbol

Map our internal YAML ``type`` string (e.g. ``"include_above"``,
``"exclude_range"``) to the canonical ``ConstraintInterval`` kind.
"""
function _interval_kind(yaml_type::AbstractString)
    sym = Symbol(yaml_type)
    sym in (:include_range, :include_above, :include_below,
            :exclude_range, :exclude_above, :exclude_below,
            :mass_estimate, :point_estimate) && return sym
    # Fallback: derive from prefix.
    startswith(yaml_type, "exclude") && return :exclude_range
    startswith(yaml_type, "include") && return :include_range
    error("Unknown interval type: $yaml_type")
end

"""
    _build_source(row) :: DataSource

Build a ``DataSource`` from a YAML row's ``bib```` dict, carrying authors,
DOI, eprint, url, title, year.
"""
function _build_source(row)
    bib = get(row, "bib", Dict{String,Any}())
    return DataSource(
        citekey = row["citekey"],
        paper   = get(bib, "title", ""),
        year    = Int(get(bib, "year", 0)),
        method  = get(row, "method", ""),
        type    = :literature,
        notes   = string("journal=", get(bib, "journal", "")),
        authors = String[String(a) for a in get(bib, "authors", String[])],
        doi     = String(get(bib, "doi", "")),
        eprint  = String(get(bib, "eprint", "")),
        bibcode = String(get(bib, "bibcode", "")),
        url     = String(get(bib, "url", "")),
    )
end

"""
    _build_intervals(row) :: Vector{ConstraintInterval}

Build one ``ConstraintInterval`` per entry in ``row["constraints"]``.
``Inf`` / ``-Inf`` are mapped from YAML ``.inf`` / ``-.inf`` (which YAML.jl
parses as ``Inf`` / ``-Inf`` already).
"""
function _build_intervals(row)
    ivs = ConstraintInterval[]
    for interval in row["constraints"]
        lo = interval["log10_min"]
        hi = interval["log10_max"]
        lo = lo === nothing ? -Inf : Float64(lo)
        hi = hi === nothing ?  Inf : Float64(hi)
        push!(ivs, ConstraintInterval(
            kind       = _interval_kind(interval["type"]),
            lower      = lo,
            upper      = hi,
            confidence = 0.95,
            label      = String(get(interval, "label", "")),
        ))
    end
    return ivs
end

"""
    row_to_constraint_row(row) :: ConstraintRow

Convert one YAML row (= one citekey) into a ``ConstraintRow`` with all its
intervals + free-text ``data`` and ``method``.
"""
function row_to_constraint_row(row)
    return ConstraintRow(
        _build_source(row),
        :log10_m_a_eV;
        intervals = _build_intervals(row),
        data      = String(get(row, "data", "")),
        method    = String(get(row, "method", "")),
    )
end

function main()
    @info "Loading YAML" path=IN_YAML
    data = YAML.load_file(IN_YAML)
    @info "Loaded citekeys" n=length(data)

    ds = DataSet{ConstraintRow}(
        "wave_cdm_constraints",
        "Observational constraints on the ultralight (wave/fuzzy) dark matter particle mass, parsed from WaveDM.jl docs/constraints/ma_constraints.tex.  Each row is one citekey with 0..N intervals plus free-text data/method labels.",
        ConstraintRow[],
    )

    for (key, row) in data
        push!(ds.data, row_to_constraint_row(row))
    end

    @info "Built DataSet" entries=length(ds)

    # -- Coverage report ----------------------------------------------------
    cat_count = Dict{String,Int}()
    for (key, row) in data
        cat = get(row, "category", "other")
        cat_count[cat] = get(cat_count, cat, 0) + 1
    end
    @info "Category coverage" cat_count

    # -- Save the DataSet --------------------------------------------------
    @info "Saving ParameterSpace YAML" path=OUT_YAML
    save_dataset(ds, OUT_YAML)

    # -- Demonstrate a query: rows from 2018-2025 ---------------------------
    filtered = ConstraintRow[r for r in ds.data if 2018 ≤ r.source.year ≤ 2025]
    sub = DataSet{ConstraintRow}(
        "wave_cdm_constraints_2018_2025",
        "Constraints published between 2018 and 2025.",
        filtered,
    )
    @info "Filtered DataSet" entries=length(sub)
    save_dataset(sub, FILTERED)

    @info "Done."
    return ds
end

main()