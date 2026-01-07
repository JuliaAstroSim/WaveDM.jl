module WaveDM

using PrecompileTools
using DocStringExtensions

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

# JuliaAstroSim encosystem
using AstroSimBase
using PhysicalParticles
using PhysicalParticles.NumericalIntegration
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

export test_MW_MOND, SPE3D_MOND
export SimulationGrid, DeviceConfig, TimeStepConfig, AstroUnitsConfig,
       GravityConfig, TidalFieldConfig, InitialConditionsConfig,
       DensityProfileConfig, MassRadiusConfig, VisualizationConfig,
       VisualizationData, ProfileFitConfig, RCFitConfig, BestFitConfig

include("core/coordinates.jl")

include("core/configs.jl")
include("core/statistics.jl")
include("core/utils.jl")
include("core/init.jl")

include("plot/plot_runtime.jl")

include("MOND.jl")

include("solver/gravity.jl")
include("solver/tidal_fields.jl")
include("solver/KDK.jl")
include("solver/best_fit_extraction.jl")

include("initial_conditions.jl")
include("io.jl")
include("simulation.jl")

include("auxiliary.jl")
include("plot/plot.jl")

include("precompile.jl")

end # module WaveDM
