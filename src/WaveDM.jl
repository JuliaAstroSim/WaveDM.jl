module WaveDM

using PrecompileTools
using DocStringExtensions

using Base.Iterators
using Random, Printf, Statistics
using Combinatorics
using Observables
using Measurements
using StatsBase

using CSV, XLSX, JLD2
using DataFrames
using StructArrays
using ColorSchemes

using Unitful, UnitfulAstro

using CUDA
using FFTW
# FFTW.set_provider!("mkl")
# FFTW.set_num_threads(Threads.nthreads())

using RollingFunctions, SpecialFunctions, Roots

using LsqFit
using Optim
using Dierckx

using Sundials

using GLMakie, CairoMakie
using UnicodePlots
using ProgressMeter

# Additional imports from temporary files
using Distributed
# using Dagger
using Strided
using GridInterpolations
using Zygote
using Tullio
using PaddedViews
using OffsetArrays
using DSP
using Dates

# Distributed-memory parallelism:
#* `PhysicalMeshes.jl` provides the adaptive `ParallelBackend`
#     abstraction (`:serial / :threads / :distributed / :gpu`).
#* `DistributedArrays.jl` is required by some algorithms (e.g. the
#     slab-decomposed 3D FFT for the distributed Poisson solver).
#* Dagger.jl is intentionally NOT loaded: see
#     JuliaParallel/Dagger.jl#649 for the DRef finalizer issue.
using DistributedArrays
using PhysicalMeshes: ParallelBackend,
                      select_backend, detect_resources,
                      is_serial, is_threads, is_distributed, is_gpu, is_parallel,
                      to_device, to_host, distribute, release!,
                      parallel_sum, parallel_maximum, parallel_minimum,
                      parallel_findmax, parallel_quantile, parallel_sumprod,
                      parallel_broadcast!, set_backend_fft_threads!

# JuliaAstroSim encosystem
using AstroSimBase
using PhysicalParticles
using PhysicalParticles.NumericalIntegration
using PhysicalParticles.BangBang
using PhysicalFDM
using PhysicalFFT
using PhysicalMeshes
using AstroIC
using AstroNbodySim
using AstroPlot
using GalacticDynamics

const C = Constant(uAstro)
astro()

const η₀ = sqrt(C.μ_0 / C.ε_0)
const ħ = C.h/2/π
const Ωₘ₀ = 0.31

export simulate_waveDM, SPE3D_waveDM
export SimulationGrid, DeviceConfig, TimeStepConfig, AstroUnitsConfig,
       GravityConfig, TidalFieldConfig, InitialConditionsConfig,
       DensityProfileConfig, MassRadiusConfig, VisualizationConfig,
       VisualizationData, ProfileFitConfig, RCFitConfig, BestFitConfig

# Beginner-facing introspection helpers — see src/ic/catalog.jl
# Each catalog gets a `list_<catalog>()` accessor (returns the AstroIC
# DataFrame) plus a `<catalog>_index(df, id)` helper (turns a 1-based user
# id into a row index, with a uniform error message).  `ufd_index` is the
# original example; the other `*_index` helpers follow the same pattern so
# adding a new catalog in `src/ic/catalog.jl` only takes a one-line
# update of the `_CATALOG_LOADERS` dict.
export list_supported_models, list_galaxies, print_catalog
export list_MW_satellites, list_SPARC_LTGs,
       list_SPARC_Xray_ETGs, list_SPARC_rotating_ETGs, list_cooke_dwarfs
export ufd_index, MW_satellite_index,
       SPARC_LTG_index, SPARC_Xray_ETG_index, SPARC_rotating_ETG_index,
       cooke_dwarf_index
# Per-catalog default `(rMin, rMax)` for RAR-style fits; see
# `src/ic/catalog.jl` for the per-catalog recipes (mirrors `SPE-SPARC.jl`).
export default_RAR_radial_window

# Re-export commonly used functions for API compatibility
export setup_grid, setup_coordinates, compute_timestep
export setup_initial_conditions, solve_vector_equation
export generate_initial_conditions
export compute_gravitational_potential, apply_kick_step!, apply_drift_step!
export setup_fft_operators, setup_absorption_boundary
export compute_profile_fit_error, compute_rc_fit_error, compute_beta_star, update_best_fit!
export compute_velocity_field_shells, compute_mass_weighted_center,
       compute_interior_mass, compute_diagnostic_snapshot,
       save_diagnostic_timeseries, detect_cuspy_windows,
       DiagnosticSnapshot
export tidal_interpolation_accuracy
export save_initial_conditions, save_evolution_results, save_property_dataframe, compute_averaged_fields
export setup_visualization, plotMOND
export need_to_interrupt, findfirstvalue
export filter_min_rho, func_dθ_dt
export vec_cartesian_to_spherical, vec_cartesian_to_cylindrical

export select_backend, detect_resources, parallel_poisson
export is_serial, is_threads, is_distributed, is_gpu, is_parallel

# Core modules (data structures and coordinates)
include("core/coordinates.jl")
include("core/configs.jl")

include("parallel_poisson.jl")

# Math utilities
include("math/statistics.jl")

# Utils modules
include("utils/kernels.jl")
include("utils/filters.jl")

# Physics modules
include("physics/MOND.jl")
include("physics/tidal_fields.jl")
include("physics/schrodinger.jl")
include("physics/poisson.jl")

# Solver modules
include("solver/gravity.jl")
include("solver/KDK.jl")
include("solver/fft_operators.jl")
include("solver/best_fit_extraction.jl")

# Initial conditions modules
include("ic/profiles.jl")
include("ic/milkyway.jl")
include("ic/catalog.jl")
include("ic/generation.jl")

# Simulation modules
include("simulation/utils.jl")
include("simulation/setup.jl")
include("simulation/loop.jl")

# IO modules
include("io/save.jl")

# Visualization modules
include("visualization/runtime.jl")
include("visualization/plots.jl")

# Main simulation files (to be split in future iterations)
include("simulation_main.jl")

include("precompile.jl")

end # module WaveDM
