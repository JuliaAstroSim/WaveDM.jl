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
using PhysicalFDM
using PhysicalFFT
using PhysicalMeshes
using AstroIC
using AstroNbodySim
using AstroPlot

export RAR, Milgrom
export RAR_inv, Milgrom_inv
export modelRAR


const C = Constant(uAstro)
astro()

const η₀ = sqrt(C.μ_0 / C.ε_0)
const ħ = C.h/2/π

include("coordinates.jl")
export vec_cartesian_to_spherical, vec_cartesian_to_cylindrical

include("statistics.jl")
include("utils.jl")
include("auxiliary.jl")

include("RAR.jl")

# Include MOND modules
include("MOND/initial_conditions.jl")
include("MOND/analysis.jl")

# Include constraints module
include("constraints.jl")

# Include high-level API module
include("SPE/init.jl")
include("SPE.jl")
include("highlevel.jl")

include("plot.jl")

include("precompile.jl")

end # module WaveDM
