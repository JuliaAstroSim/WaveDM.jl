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
using Dierckx
using DifferentialEquations
using ProgressMeter

using Sundials

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

export load_SPARC_LTGs_RC, load_SPARC_LTGs_data, load_li2018_SPARC
export load_SPARC_Xray_ETGs_data, load_SPARC_rotating_ETGs_data, load_SPARC_rotating_ETGs_RC, load_SPARC_rotating_ETGs_rotmod


const C = Constant(uAstro)
astro()

const η₀ = sqrt(C.μ_0 / C.ε_0)
const ħ = C.h/2/π


include("data.jl")
include("statistics.jl")

include("RAR.jl")
include("SPE.jl")

include("plot.jl")

include("precompile.jl")

end # module WaveDM
