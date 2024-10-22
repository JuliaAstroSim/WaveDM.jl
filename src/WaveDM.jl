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

export load_SPARC_RC, load_SPARC_data, load_li2018_SPARC


include("data.jl")

include("RAR.jl")
include("SPE.jl")


include("precompile.jl")

end # module WaveDM
