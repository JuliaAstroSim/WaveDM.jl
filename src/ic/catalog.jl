"""
    list_supported_models()

Return a `Vector{NamedTuple}` describing every `model = :...` keyword accepted
by [`simulate_waveDM`](@ref). Each entry has the fields

| field        | meaning                                                 |
| ------------ | ------------------------------------------------------- |
| `model`      | the symbol you pass to the `model` keyword              |
| `baryon_mode`| `Vector{Symbol}` listing the allowed `baryon_mode`s    |
| `description`| one-line description of what the model represents       |

This is the **first thing beginners should run** when they want to know which
galaxies / halos are built into WaveDM.jl.
"""
function list_supported_models()
    return [
        (
            model = :MW,
            baryon_mode = [:ignored, :mesh, :particles_static, :particles_dynamic],
            description = "Milky Way analogue (Zhu 2023 mass model: bulge + thin/thick stellar disc + HI/HII gas + NFW DM halo).",
        ),

        (
            model = :dwarf_UFDs,
            baryon_mode = [:ignored],
            description = "Ultra-faint dwarf galaxies (Crater II, Draco, Ursa Minor, …). Data from Hayashi et al. 2023 via AstroIC.load_data_UFDs(); `Galaxy_id` is the 1-based row index into that table (Crater II is row 6).",
        ),

        (
            model = :dwarf,
            baryon_mode = [:ignored, :particles_static, :particles_dynamic],
            description = "Generic dwarf galaxy with a gNFW halo and an optional stellar + gaseous disc.",
        ),

        (
            model = :dwarf_NFW,
            baryon_mode = [:ignored],
            description = "Generic dwarf galaxy with a classical NFW halo (β = 1).",
        ),

        (
            model = :dwarf_Zhao,
            baryon_mode = [:ignored],
            description = "Generic dwarf galaxy with a generalized Zhao (1996) profile (free α, β, γ).",
        ),

        (
            model = :Elliptical,
            baryon_mode = [:mesh, :ignored],
            description = "Early-type galaxy with a gNFW halo and a Jaffe-model stellar component on the mesh.",
        ),

        (
            model = :cluster_NFW,
            baryon_mode = [:mesh, :ignored],
            description = "Galaxy cluster with a gNFW halo and a β-model ICM (intracluster medium) on the mesh.",
        ),

        (
            model = :cluster_Burkert,
            baryon_mode = [:mesh, :ignored],
            description = "Galaxy cluster with a Burkert (1992) core halo and a β-model ICM.",
        ),

        (
            model = :SPARC_LTGs,
            baryon_mode = [:ignored, :mesh, :particles_static, :particles_dynamic],
            description = "SPARC late-type galaxies (planned).",
        ),

        (
            model = :SPARC_Xray_ETGs,
            baryon_mode = [:ignored, :mesh],
            description = "SPARC X-ray early-type galaxies (planned).",
        ),

        (
            model = :SPARC_rotating_ETGs,
            baryon_mode = [:ignored, :mesh],
            description = "SPARC rotating early-type galaxies (planned).",
        ),
    ]
end

# ===========================================================================
# Catalog dispatch tables
# ---------------------------------------------------------------------------
# Every WaveDM "named-object" model (`:dwarf_UFDs`, `:MW`, `:SPARC_*`,
# `:dwarf_massive`) pulls its object list from an AstroIC.jl loader:
#
#   * `AstroIC.load_data_UFDs()`             — Hayashi 2023 UFDs
#   * `AstroIC.load_data_MW_satellites()`    — Battaglia 2022 MW satellites
#   * `AstroIC.load_SPARC_LTGs_data()`       — SPARC late-type galaxies
#   * `AstroIC.load_SPARC_Xray_ETGs_data()`  — SPARC X-ray ETGs
#   * `AstroIC.load_SPARC_rotating_ETGs_data()` — SPARC rotating ETGs
#   * `AstroIC.load_massive_dwarf_CO_RC`     — Cooke 2022 massive dwarfs
#
# To stay flexible as new catalogs are added, we expose one
# `list_<catalog>()` function per catalog that returns the underlying
# DataFrame, plus one `<catalog>_index(df, id)` helper per catalog that
# translates a user-facing 1-based id into the same DataFrame (with a
# descriptive error if `id` is out of range).  WaveDM never owns a hardcoded
# table; everything comes from AstroIC.
# ===========================================================================

# All known catalog loaders, keyed by their model symbol.  Used by
# `list_galaxies` (and the future `:MW`, `:SPARC_*`, … dispatchers) so
# adding a new catalog is a one-line change here.
const _CATALOG_LOADERS = Dict{Symbol,Function}(
    :dwarf_UFDs           => AstroIC.load_data_UFDs,
    :MW_satellites        => AstroIC.load_data_MW_satellites,
    :SPARC_LTGs           => AstroIC.load_SPARC_LTGs_data,
    :SPARC_Xray_ETGs      => AstroIC.load_SPARC_Xray_ETGs_data,
    :SPARC_rotating_ETGs  => AstroIC.load_SPARC_rotating_ETGs_data,
)

# The Cooke 2022 catalogs are per-galaxy files, not a single DataFrame, so
# they get their own accessor (see `list_cooke_dwarfs` below).
const _COOKE_DWARFS = ("NGC1035", "NGC4310", "NGC4451", "NGC4701",
                       "NGC5692", "NGC6106")

"""
    _out_of_range_error(label, model, id, n_rows, df, name_col)

Internal helper used by every `*_index` function below to throw a uniform
"id out of range" error.  `label` is e.g. `"Galaxy_id"`, `model` is e.g.
`:dwarf_UFDs`, `id` is the user-supplied integer, `n_rows` is `nrow(df)`,
`df` is the underlying DataFrame, and `name_col` is the symbol of the
column whose values are listed next to each id (e.g. `:Galaxy`).
"""
function _out_of_range_error(label::AbstractString, model::Symbol,
                              id::Integer, n_rows::Integer,
                              df, name_col::Symbol)
    supported = string.(1:n_rows) .* " (" .* string.(df[!, name_col]) .* ")"
    loader = get(_CATALOG_LOADERS, model, nothing)
    loader_msg = if loader === nothing
        "the catalog table"
    else
        "`AstroIC.$(string(nameof(loader)))`"
    end
    msg = """
    `$(label) = $(id)` is out of range for `model = :$(model)`.

    `$(label)` is the 1-based row index into the table loaded by
    $(loader_msg).  Valid values:
        $(join(supported, ", "))

    See `WaveDM.list_galaxies()` (or the catalog-specific helper) for the
    full table.

    In the meantime, you can fall back to a generic profile
    (`model = :dwarf_NFW` / `:dwarf_Zhao`) and pass `halo_ρ0` / `halo_r0`
    / `halo_β` via `config_profile = ...`.
    """
    error(msg)
end

# ===========================================================================
# UFD catalog — Hayashi et al. 2023 (Tables 1 + A.1)
# ===========================================================================

"""
    list_galaxies(; model = :dwarf_UFDs)

Return the `DataFrame` of available objects for `model`.  The columns come
straight from `AstroIC.load_data_UFDs()`:

| column      | unit     | meaning                                    |
| ----------- | -------- | ------------------------------------------ |
| `Galaxy`    | —        | common name (e.g. "Crater_2")              |
| `N_samples` | —        | number of stars in the Hayashi sample      |
| `L`         | L☉       | V-band luminosity (10^log10L)              |
| `b`         | pc       | projected half-light radius (10^log10b)    |
| `rho0`      | Msun/pc³ | central DM density (10^log10rho0)          |
| `alpha`     | —        | Zhao-profile α fit                         |
| `beta`      | —        | Zhao-profile β fit (inner slope)           |
| `gamma`     | —        | Zhao-profile γ fit (outer slope)           |
| `Q`,`Q_d`,`Q_u` | —    | anisotropy parameter + errors               |

Every column also has a `_u` / `_d` companion giving the upper / lower
uncertainty on the logarithm.

If you don't know which `Galaxy_id` to use, run this function first and
pick the row whose `Galaxy` matches the dwarf you want.

## Example

```julia
julia> using WaveDM

julia> df = list_galaxies()
```
"""
function list_galaxies(; model::Symbol = :dwarf_UFDs)
    loader = get(_CATALOG_LOADERS, model, nothing)
    loader === nothing && throw(ArgumentError(
        "`list_galaxies(; model = :$model)` is not wired up in this release. " *
        "For other models pass parameters through `config_profile = ...` instead. " *
        "Run `list_supported_models()` to see what is currently implemented."
    ))
    return loader()
end

"""
    ufd_index(df_UFDs, Galaxy_id) -> row_index

Translate the user-facing `Galaxy_id` (1-based row index into the Hayashi
2023 table loaded by [`AstroIC.load_data_UFDs`](@ref)) into a row index
that downstream code (`get_dwarf_UFDs_params`, etc.) can use directly.

Throws a uniform "out of range" error (powered by
[`_out_of_range_error`](@ref)) if `Galaxy_id` is outside `[1, nrow(df_UFDs)]`.
"""
function ufd_index(df_UFDs, Galaxy_id::Integer)
    n = nrow(df_UFDs)
    if Galaxy_id < 1 || Galaxy_id > n
        _out_of_range_error("Galaxy_id", :dwarf_UFDs, Int(Galaxy_id), n, df_UFDs, :Galaxy)
    end
    return Int(Galaxy_id)
end

# ===========================================================================
# MW satellites — Battaglia et al. 2022 (Gaia EDR3 systemic motions)
# ===========================================================================

"""
    list_MW_satellites()

Convenience accessor for `AstroIC.load_data_MW_satellites()`.  Returns the
DataFrame of Milky-Way satellite galaxies (Sagittarius, LMC, SMC, Draco,
…); see [Battaglia et al. 2022](https://www.aanda.org/) for the source.
"""
list_MW_satellites() = AstroIC.load_data_MW_satellites()

"""
    MW_satellite_index(df_MW_satellites, Galaxy_id) -> row_index

Translate a user-facing `Galaxy_id` (1-based row index into the Battaglia
2022 table loaded by [`AstroIC.load_data_MW_satellites`](@ref)) into a row
index for downstream code.

Throws a uniform "out of range" error if `Galaxy_id` is outside
`[1, nrow(df_MW_satellites)]`.
"""
function MW_satellite_index(df_MW_satellites, Galaxy_id::Integer)
    n = nrow(df_MW_satellites)
    if Galaxy_id < 1 || Galaxy_id > n
        _out_of_range_error("Galaxy_id", :MW_satellites, Int(Galaxy_id), n,
                            df_MW_satellites, :Galaxy)
    end
    return Int(Galaxy_id)
end

# ===========================================================================
# SPARC late-type galaxies — Lelli et al. 2016c
# ===========================================================================

"""
    list_SPARC_LTGs()

Convenience accessor for `AstroIC.load_SPARC_LTGs_data()`.  Returns the
filtered DataFrame of late-type galaxies from the SPARC sample (Lelli
et al. 2016c); see `src/data/SPARC.jl` for the applied cuts.
"""
list_SPARC_LTGs() = AstroIC.load_SPARC_LTGs_data()

"""
    SPARC_LTG_index(df_SPARC_LTGs, Galaxy_id) -> row_index

Translate a user-facing `Galaxy_id` (1-based row index into the SPARC
late-type galaxy table) into a row index for downstream code.  Throws a
uniform "out of range" error if `Galaxy_id` is outside
`[1, nrow(df_SPARC_LTGs)]`.
"""
function SPARC_LTG_index(df_SPARC_LTGs, Galaxy_id::Integer)
    n = nrow(df_SPARC_LTGs)
    if Galaxy_id < 1 || Galaxy_id > n
        _out_of_range_error("Galaxy_id", :SPARC_LTGs, Int(Galaxy_id), n,
                            df_SPARC_LTGs, :Galaxy)
    end
    return Int(Galaxy_id)
end

# ===========================================================================
# SPARC X-ray ETGs — Lelli et al. 2017
# ===========================================================================

"""
    list_SPARC_Xray_ETGs()

Convenience accessor for `AstroIC.load_SPARC_Xray_ETGs_data()`.  Returns
the X-ray early-type galaxy sample from Lelli et al. 2017.
"""
list_SPARC_Xray_ETGs() = AstroIC.load_SPARC_Xray_ETGs_data()

"""
    SPARC_Xray_ETG_index(df, Galaxy_id) -> row_index

Translate a user-facing `Galaxy_id` (1-based row index into the SPARC
X-ray ETG table) into a row index for downstream code.  Throws a uniform
"out of range" error if `Galaxy_id` is outside `[1, nrow(df)]`.
"""
function SPARC_Xray_ETG_index(df, Galaxy_id::Integer)
    n = nrow(df)
    if Galaxy_id < 1 || Galaxy_id > n
        _out_of_range_error("Galaxy_id", :SPARC_Xray_ETGs, Int(Galaxy_id), n,
                            df, :Galaxy)
    end
    return Int(Galaxy_id)
end

# ===========================================================================
# SPARC rotating ETGs — Lelli et al. 2017
# ===========================================================================

"""
    list_SPARC_rotating_ETGs()

Convenience accessor for `AstroIC.load_SPARC_rotating_ETGs_data()`.  Returns
the rotating early-type galaxy sample from Lelli et al. 2017.
"""
list_SPARC_rotating_ETGs() = AstroIC.load_SPARC_rotating_ETGs_data()

"""
    SPARC_rotating_ETG_index(df, Galaxy_id) -> row_index

Translate a user-facing `Galaxy_id` (1-based row index into the SPARC
rotating ETG table) into a row index for downstream code.  Throws a uniform
"out of range" error if `Galaxy_id` is outside `[1, nrow(df)]`.
"""
function SPARC_rotating_ETG_index(df, Galaxy_id::Integer)
    n = nrow(df)
    if Galaxy_id < 1 || Galaxy_id > n
        _out_of_range_error("Galaxy_id", :SPARC_rotating_ETGs, Int(Galaxy_id), n,
                            df, :Galaxy)
    end
    return Int(Galaxy_id)
end

# ===========================================================================
# Cooke 2022 massive dwarfs — per-galaxy RC CSVs (NGC1035, NGC4310, …)
# ===========================================================================

"""
    list_cooke_dwarfs()

Return the tuple of galaxy names available in the Cooke 2022 massive-dwarf
rotation-curve sample.  Each name maps to three CSV files
(`Cooke2022_RC/<name>.csv`, `Cooke2022_RC_DM/<name>.csv`,
`Cooke2022_RC_stellar/<name>.csv`); use [`load_massive_dwarf_CO_RC`](@ref),
[`load_massive_dwarf_DM_RC`](@ref), and [`load_massive_dwarf_Baryon_RC`](@ref)
to access them.

This is the Cooke-side equivalent of `list_galaxies()` — instead of a
single DataFrame we hand back a tuple because each galaxy is its own
file.
"""
list_cooke_dwarfs() = _COOKE_DWARFS

"""
    cooke_dwarf_index(name::AbstractString) -> Galaxy_id

Translate a galaxy name into the 1-based index used elsewhere
(`plot_RC_RAR(..., Galaxy_id, ...)` etc.).  Throws a uniform "out of range"
error if `name` is not in the Cooke 2022 sample.
"""
function cooke_dwarf_index(name::AbstractString)
    id = findfirst(==(name), _COOKE_DWARFS)
    id === nothing && error(
        "`$(name)` is not in the Cooke 2022 massive-dwarf sample. " *
        "Available names: $(join(_COOKE_DWARFS, ", "))."
    )
    return id
end

# ===========================================================================
# Per-catalog RAR radial-window defaults
# ---------------------------------------------------------------------------
# The per-snapshot RAR a_0 diagnostic in `SPE3D_waveDM` needs a radial
# window `[rMin, rMax]`.  Crater II (~1 kpc) and the MW (~30 kpc) need
# *very* different windows, so it pays to look up the catalog-appropriate
# defaults from the same place where the catalog lives.  These helpers
# mirror the choices used in `SPE-SPARC.jl`:
#
#   * SPARC X-ray ETGs           — `rMin = 0.1 × Reff`,  `rMax = 4 × Reff`
#   * SPARC rotating ETGs        — `rMax = 4 × max(Reff, Rexp)`
#   * SPARC LTGs                 — `rMax = 2 × max(rs, RHI)`
#   * UFDs (Hayashi 2023)        — `rMin = b / 2`,       `rMax = 4 × b`
#                                  where `b` is the projected half-light
#                                  radius (the natural "Reff" analog for
#                                  dwarf galaxies)
#   * MW                         — fall back to the hardcoded
#                                  `minR = 8 kpc` / `maxR = 15 kpc` of
#                                  the SPE driver.
#
# The function takes a `catalog::Symbol` (`:dwarf_UFDs`, `:SPARC_LTGs`,
# `:SPARC_Xray_ETGs`, `:SPARC_rotating_ETGs`, or `:MW`) and an optional
# `Galaxy_id` (1-based row index into the matching DataFrame), and returns
# a `(minR_kpc, maxR_kpc)::Tuple{Float64,Float64}`.
# ===========================================================================

"""
    default_RAR_radial_window(catalog::Symbol;
                              Galaxy_id::Integer = 1,
                              minR_kpc::Real = 0.0,
                              maxR_kpc::Real = Inf) -> (minR_kpc, maxR_kpc)

Return a per-catalog default `(rMin, rMax)` window (both in kpc) suitable
for fitting the Radial Acceleration Relation (RAR) for the galaxy at
`Galaxy_id` in the given `catalog`.

`minR_kpc` / `maxR_kpc` act as floor / ceiling overrides — if the catalog
default would be smaller than `minR_kpc`, that is returned; if larger than
`maxR_kpc`, that is returned.  This matches the `max(rMin, rMin_cell_i, 0.1 Reff)`
clamp pattern used in `SPE-SPARC.jl`.

| `catalog`            | source                                              |
| -------------------- | --------------------------------------------------- |
| `:dwarf_UFDs`        | Hayashi 2023 — `rMin = b/2`, `rMax = 4 × b` (kpc)   |
| `:SPARC_LTGs`        | Lelli 2016c — `rMax = 2 × max(rs, RHI)` (kpc)       |
| `:SPARC_Xray_ETGs`   | Lelli 2017 — `rMin = 0.1 × Reff`, `rMax = 4 × Reff` |
| `:SPARC_rotating_ETGs` | Lelli 2017 — `rMax = 4 × max(Reff, Rexp)`           |
| `:MW`                | falls back to `(minR_kpc, maxR_kpc)` (the SPE default) |

For other catalogs (`:dwarf_massive`, etc.) the helper falls back to
`(minR_kpc, maxR_kpc)` — users can override by passing their own.
"""
function default_RAR_radial_window(catalog::Symbol;
                                    Galaxy_id::Integer = 1,
                                    minR_kpc::Real = 0.0,
                                    maxR_kpc::Real = Inf)
    base = (Float64(minR_kpc), isinf(maxR_kpc) ? 0.0 : Float64(maxR_kpc))
    catalog === :MW && return base

    if catalog === :dwarf_UFDs
        df = AstroIC.load_data_UFDs()
        1 ≤ Galaxy_id ≤ nrow(df) || throw(ArgumentError(
            "`Galaxy_id = $Galaxy_id` is out of range for `:dwarf_UFDs` (1..$(nrow(df)))."
        ))
        # Hayashi 2023 `b` is the projected half-light radius in pc.
        # Use `b/2` as the inner cutoff (consistent with `rMin = 0.1 Reff`
        # for Reff ≈ 5 b) and `4 b` as the outer cutoff (analogous to
        # `rMax = 4 Reff` in the SPARC ETG routines).
        b_kpc = ustrip(u"kpc", df.b[Galaxy_id])
        rmin_cat = b_kpc / 2
        rmax_cat = 4 * b_kpc
    elseif catalog === :SPARC_LTGs
        df = AstroIC.load_SPARC_LTGs_data()
        1 ≤ Galaxy_id ≤ nrow(df) || throw(ArgumentError(
            "`Galaxy_id = $Galaxy_id` is out of range for `:SPARC_LTGs` (1..$(nrow(df)))."
        ))
        # SPARC LTG column names are `rs` (scale radius, kpc) and `RHI`
        # (HI disk radius, kpc).  Both are kpc in the loader.
        rs_kpc  = Float64(df.rs[Galaxy_id])
        RHI_kpc = Float64(df.RHI[Galaxy_id])
        rmax_cat = 2 * max(rs_kpc, RHI_kpc)
        rmin_cat = 0.0  # SPARC LTGs do not impose a fixed rMin in SPE-SPARC.jl
    elseif catalog === :SPARC_Xray_ETGs
        df = AstroIC.load_SPARC_Xray_ETGs_data()
        1 ≤ Galaxy_id ≤ nrow(df) || throw(ArgumentError(
            "`Galaxy_id = $Galaxy_id` is out of range for `:SPARC_Xray_ETGs` (1..$(nrow(df)))."
        ))
        Reff_kpc = Float64(df.Reff[Galaxy_id])
        rmin_cat = 0.1 * Reff_kpc
        rmax_cat = 4.0  * Reff_kpc
    elseif catalog === :SPARC_rotating_ETGs
        df = AstroIC.load_SPARC_rotating_ETGs_data()
        1 ≤ Galaxy_id ≤ nrow(df) || throw(ArgumentError(
            "`Galaxy_id = $Galaxy_id` is out of range for `:SPARC_rotating_ETGs` (1..$(nrow(df)))."
        ))
        Reff_kpc = Float64(df.Reff[Galaxy_id])
        Rexp_kpc = Float64(df.Rexp[Galaxy_id])
        rmax_cat = 4 * max(Reff_kpc, Rexp_kpc)
        rmin_cat = 0.0  # SPE-SPARC.jl uses no explicit rMin for rotating ETGs
    else
        # Unknown / unsupported catalog: fall back to the user-supplied
        # floor / ceiling unchanged.
        return base
    end

    # Floor / ceiling overrides
    rmin = max(rmin_cat, Float64(minR_kpc))
    rmax = if isinf(maxR_kpc)
        rmax_cat
    else
        min(rmax_cat, Float64(maxR_kpc))
    end
    return (rmin, rmax)
end

"""
    _resolve_rar_window(minR_kpc, maxR_kpc, catalog::Symbol; Galaxy_id = 1,
                       floor_minR_kpc = 0.0, cap_maxR_kpc = Inf)

Resolve a user-supplied `(minR_kpc, maxR_kpc)` pair into concrete kpc
floats (all returned values are plain `Float64`).  Each of `minR_kpc` /
`maxR_kpc` may be either a numeric `Real` (in kpc) or the symbol `:auto` —
in the latter case the per-catalog default from
[`default_RAR_radial_window`](@ref) is used, with the *other* override
acting as the floor / ceiling.

`floor_minR_kpc` and `cap_maxR_kpc` follow the same convention: any
`Quantity{<:Any, 𝐋, ...}` is stripped to kpc via `ustrip(u"kpc", x)`.
Callers in `SPE3D_waveDM` strip units before calling (using the same
`ustrip(u"kpc", …)` pattern as `uL = ustrip(u"kpc", length_astro)`)
but the helper stays friendly to bare-`Real` and `Quantity` inputs.

This is the helper used by [`SPE3D_waveDM`](@ref)'s `flag_RAR_snapshot`
path so users can write `RAR_snapshot_minR_kpc = :auto, RAR_snapshot_maxR_kpc = :auto`
and get sensible defaults for any catalog / galaxy combination.
"""
function _resolve_rar_window(minR_kpc, maxR_kpc, catalog::Symbol;
                              Galaxy_id::Integer = 1,
                              floor_minR_kpc = 0.0,
                              cap_maxR_kpc = Inf)
    # Tolerate either a bare Real or a Quantity[<length unit>] — strip to kpc.
    _kpc(x::Real)    = Float64(x)
    _kpc(x::Quantity) = Float64(ustrip(u"kpc", x))
    _kpc(::Nothing)  = 0.0

    auto = :auto
    is_auto_min = minR_kpc === auto || (minR_kpc isa Real && minR_kpc < 0)
    is_auto_max = maxR_kpc === auto || (maxR_kpc isa Real && maxR_kpc < 0)

    if is_auto_min || is_auto_max
        cat_min, cat_max = default_RAR_radial_window(catalog;
                                                      Galaxy_id = Galaxy_id,
                                                      minR_kpc = 0.0,
                                                      maxR_kpc = Inf)
        rmin = is_auto_min ? cat_min : _kpc(minR_kpc)
        rmax = is_auto_max ? cat_max : _kpc(maxR_kpc)
    else
        rmin = _kpc(minR_kpc)
        rmax = _kpc(maxR_kpc)
    end
    # Apply the universal floor / ceiling from the caller's kwarg defaults
    rmin = max(rmin, _kpc(floor_minR_kpc))
    rmax = if isinf(_kpc(cap_maxR_kpc))
        rmax
    else
        min(rmax, _kpc(cap_maxR_kpc))
    end
    return (rmin, rmax)
end

# ===========================================================================
# Human-readable catalog summary
# ===========================================================================

"""
    print_catalog()

Print a human-readable summary of every supported model and every available
catalog object, formatted for the REPL.  Useful as a one-liner cheatsheet:

```julia
julia> using WaveDM; print_catalog()
```
"""
function print_catalog()
    println("=" ^ 70)
    println("WaveDM.jl — supported models")
    println("=" ^ 70)
    for m in list_supported_models()
        bm = join(string.(m.baryon_mode), ", ")
        println()
        println("  model       : $(m.model)")
        println("  baryon_mode : $(bm)")
        println("  description : $(m.description)")
    end

    # ---------- UFDs (Hayashi 2023) ----------------------------------------
    println()
    println("=" ^ 70)
    println("Available ultra-faint dwarf galaxies (model = :dwarf_UFDs)")
    println("=" ^ 70)
    println("  Source: Hayashi et al. 2023 via AstroIC.load_data_UFDs().")
    println("  `Galaxy_id` is the 1-based row index into that table.")
    println()
    df = list_galaxies()
    @printf("  %-3s  %-22s  %-12s  %-10s  %-10s\n",
            "id", "name", "L [Lsun]", "b [pc]", "rho0 [Msun/pc^3]")
    for i in 1:nrow(df)
        @printf("  %-3d  %-22s  %-12.3e  %-10.1f  %-12.3e\n",
                i, df.Galaxy[i],
                ustrip(df.L[i]),
                ustrip(u"pc", df.b[i]),
                ustrip(u"Msun/pc^3", df.rho0[i]))
    end

    # ---------- MW satellites (Battaglia 2022) ----------------------------
    println()
    println("=" ^ 70)
    println("Milky-Way satellites (model = :MW — Battaglia et al. 2022)")
    println("=" ^ 70)
    df_mw = list_MW_satellites()
    @printf("  %-3s  %-15s  %-10s  %-12s\n",
            "id", "name", "d [kpc]", "v_los [km/s]")
    for i in 1:nrow(df_mw)
        @printf("  %-3d  %-15s  %-10.1f  %-+12.2f\n",
                i, df_mw.Galaxy[i],
                ustrip(u"kpc", df_mw.Distance[i]),
                ustrip(u"km/s", df_mw.v_los[i]))
    end

    # ---------- SPARC LTGs (Lelli 2016c) ----------------------------------
    println()
    println("=" ^ 70)
    println("SPARC late-type galaxies (planned: model = :SPARC_LTGs)")
    println("=" ^ 70)
    df_lt = list_SPARC_LTGs()
    @printf("  %-3s  %-12s  %-8s  %-8s\n",
            "id", "name", "D [Mpc]", "Vflat")
    for i in 1:nrow(df_lt)
        @printf("  %-3d  %-12s  %-8.2f  %-8.1f\n",
                i, df_lt.Galaxy[i], df_lt.D[i], df_lt.Vflat[i])
    end

    # ---------- SPARC X-ray ETGs (Lelli 2017) -----------------------------
    println()
    println("=" ^ 70)
    println("SPARC X-ray early-type galaxies (planned: model = :SPARC_Xray_ETGs)")
    println("=" ^ 70)
    df_x = list_SPARC_Xray_ETGs()
    @printf("  %-3s  %-12s  %-8s  %-8s\n",
            "id", "name", "D [Mpc]", "L")
    for i in 1:nrow(df_x)
        @printf("  %-3d  %-12s  %-8.2f  %-8.2e\n",
                i, df_x.Galaxy[i], df_x.D[i], df_x.L[i])
    end

    # ---------- SPARC rotating ETGs (Lelli 2017) --------------------------
    println()
    println("=" ^ 70)
    println("SPARC rotating early-type galaxies (planned: model = :SPARC_rotating_ETGs)")
    println("=" ^ 70)
    df_r = list_SPARC_rotating_ETGs()
    @printf("  %-3s  %-12s  %-8s  %-8s\n",
            "id", "name", "Dist", "Inc")
    for i in 1:nrow(df_r)
        @printf("  %-3d  %-12s  %-8.2f  %-8.1f\n",
                i, df_r.Galaxy[i], df_r.Dist[i], df_r.Inc[i])
    end

    # ---------- Cooke 2022 massive dwarfs ---------------------------------
    println()
    println("=" ^ 70)
    println("Cooke 2022 massive dwarfs (model = :dwarf_massive)")
    println("=" ^ 70)
    cooke = list_cooke_dwarfs()
    for (i, name) in enumerate(cooke)
        println("  $i  $name")
    end

    println()
    return nothing
end

# ===========================================================================
# REPL-friendly aliases
# ===========================================================================
#
#   julia> WaveDM.models
#   julia> WaveDM.galaxies
#
const models    = list_supported_models
const galaxies  = list_galaxies

# ---------------------------------------------------------------------------
# Legacy alias kept for backwards compatibility with `src/ic/generation.jl`
# and downstream scripts (SPE-WaveDM examples, …) that used the old
# DWARF_UFD_TABLE-only error helper.  All new code should use
# `ufd_index(...)` which throws via `_out_of_range_error(...)`.
# ---------------------------------------------------------------------------
function _unsupported_ufd_error(id::Int)
    df = AstroIC.load_data_UFDs()
    _out_of_range_error("Galaxy_id", :dwarf_UFDs, id, nrow(df), df, :Galaxy)
end