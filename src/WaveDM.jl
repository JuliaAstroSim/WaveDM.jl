module WaveDM

using PrecompileTools
using DocStringExtensions

using Random, Printf, Statistics
using Combinatorics
using Observables
using Measurements

using CSV, XLSX, JLD2
using DataFrames, StructArrays
using ColorSchemes

using Unitful, UnitfulAstro

using CUDA
using FFTW

using RollingFunctions, SpecialFunctions, Roots

using LsqFit
using Optim
using Dierckx
# using DifferentialEquations
using ProgressMeter

using Sundials

using GLMakie

# Additional imports from temporary files
using Distributed
using Dagger
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

const C = Constant(uAstro)
astro()

const η₀ = sqrt(C.μ_0 / C.ε_0)
const ħ = C.h/2/π

include("core/coordinates.jl")

include("core/statistics.jl")
include("core/utils.jl")
include("core/init.jl")

include("initial_conditions.jl")
include("plot/plot_runtime.jl")

include("MOND.jl")

include("solver/gravity.jl")
include("solver/tidal_fields.jl")
include("solver/KDK.jl")
include("solver/best_fit_extraction.jl")

include("initial_conditions.jl")
include("simulation.jl")

include("auxiliary.jl")
include("plot/plot.jl")

include("precompile.jl")

end # module WaveDM
