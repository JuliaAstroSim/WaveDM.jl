# Quickstart

This page is for **first-time users**. If you already know what `simulate_waveDM` is and you want to see all the knobs, jump straight to [Examples](@ref) or [Algorithms](@ref).

## 0 · Sanity-check the install (30 s)

```julia
julia> using WaveDM
julia> WaveDM.print_catalog()
```

If you see a table of supported models and ultra-faint dwarf galaxies, the install worked. If you see a `Package WaveDM not found`, go back to [Installation](@ref).

## 1 · What's in WaveDM.jl — pick your model

Before you write any code, run:

```julia
julia> list_supported_models()
```

You'll get back a list of every `model = :...` keyword accepted by `simulate_waveDM`, with three useful columns:

| `kind`       | meaning                                                    |
| ------------ | ---------------------------------------------------------- |
| `"halo"`     | big dark-matter halo (Milky-Way analogue, SPARC later)     |
| `"dwarf"`    | dwarf galaxy halo (Crater II, Draco, …)                    |
| `"cluster"`  | galaxy cluster                                             |
| `"elliptical"` | early-type galaxy                                        |

Each entry also lists the allowed `baryon_mode` values. If you don't know what `baryon_mode` is, see [Concepts below](#4-baryon_modes-in-one-line).

## 2 · Pick a galaxy (or stay generic)

If you want a *named* galaxy (Crater II, Draco, Fornax, …), use the ultra-faint-dwarf preset:

```julia
julia> list_galaxies()        # also accessible as WaveDM.galaxies()
```

This prints every `Galaxy_id` you can plug into `simulate_waveDM(; model = :dwarf_UFDs, Galaxy_id = N)`.

If you want a *custom* halo (any density profile, any parameters), skip the preset and use one of the generic models:

| If you want to …                          | Use this model        |
| ----------------------------------------- | --------------------- |
| try a Milky-Way-sized halo                | `:MW`                 |
| try a dwarf galaxy with a generic gNFW    | `:dwarf`              |
| try a classical NFW dwarf (β = 1)         | `:dwarf_NFW`          |
| try a Zhao (1996) profile (free α, β, γ)  | `:dwarf_Zhao`         |
| try an early-type galaxy                  | `:Elliptical`         |
| try a galaxy cluster                      | `:cluster_NFW` or `:cluster_Burkert` |

## 3 · The shortest script that actually runs

Copy-paste this into a file `crater2.jl` and run `julia --project crater2.jl`:

```julia
using WaveDM
using Unitful

simulate_waveDM(;
    model            = :dwarf_UFDs,
    Galaxy_id        = 6,                    # Crater II
    V                = (x, y, z, ψ) -> 0.0,  # no external potential
    Nx               = 128,                  # try 384 for production
    Xmax             = 20u"kpc",             # half-side length
    Tmax             = 1.0u"Gyr",            # simulation duration
    autoset_timestep = true,
    gpu              = false,                # set true on a GPU machine
    Realtime         = false,                # set true if you have a monitor
    title            = "MyFirstRun",
    outputdir        = "results/MyFirstRun",
)
```

What this does: builds a wave-CDM Crater II halo (gNFW, β = 1, $r_s = 3$ kpc, Hayashi 2023 defaults), evolves it for 1 Gyr on a 128³ Cartesian grid in a 40 kpc box, writes a snapshot series and a CSV of total mass / mass-fraction radii / virial energies to `results/MyFirstRun/`.

## 4 · `baryon_mode` in one line

`baryon_mode` tells WaveDM.jl how to handle the *non-dark-matter* mass in your system:

| `baryon_mode`         | What it does                                                  | When to use                                |
| --------------------- | ------------------------------------------------------------- | ------------------------------------------ |
| `:ignored`            | No baryons at all.                                            | DM-only halos (Crater II, isolated dwarfs).|
| `:mesh`               | Baryons placed on the same grid as the wave field.            | Bulges, ICM, simple extended baryons.      |
| `:particles_static`   | Baryons as N-body particles, **frozen** during the SPE loop.  | Realistic MW-like host potentials.         |
| `:particles_dynamic`  | Baryons as N-body particles, **co-evolving** with the wave field. | Full coupling, but expensive (~10× more). |

Use `:ignored` until you have a reason not to. Each model only accepts a subset — `list_supported_models()` tells you which.

## 5 · Common beginner mistakes

- **"Unknown model: :dwarf_UFDs"** → Make sure you're on the latest release of WaveDM.jl. If the model name is one of `:SPARC_LTGs / :SPARC_Xray_ETGs / :SPARC_rotating_ETGs`, those are planned but not implemented yet.
- **"Galaxy_id = 5 is out of range"** → `Galaxy_id` is the 1-based row index into the Hayashi 2023 UFD table (Crater II is row 6).  Run `list_galaxies()` to see valid values.  For other dwarfs you can use `model = :dwarf_NFW` / `:dwarf_Zhao` and pass `config_profile = ...`, or use `model = :dwarf_UFDs` with `Galaxy_id = 6` for Crater II.
- **"Out of memory" at `Nx = 384`** → That's a 384³ complex array ≈ 4 GB on CPU. Drop to `Nx = 192`, or move to a GPU.
- **Nothing happens for a long time at startup** → `simulate_waveDM` is building initial conditions (gNFW sampling + FFT-Poisson solve). Watch the `@info` lines; if you don't see them, your call didn't dispatch — check `model = :...` spelling.
- **Baryon particles complain about `Np`** → The model needs you to pass `Np` for particle counts. `list_supported_models()` shows `needs_Np = yes/no`.

## 6 · Where to go next

- [Examples](@ref) — 10 worked scripts, from a 30-second smoke test to Crater II with the MW + LMC tidal field.
- [Algorithms](@ref) — the math behind the SPE, the KDK integrator, the resolution criteria.
- [API / configs](@ref) — every keyword argument of `simulate_waveDM`.