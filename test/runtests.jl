using Test
using FFTW
using Unitful, UnitfulAstro

using AstroSimBase
using PhysicalParticles
using AstroNbodySim

using WaveDM

using CairoMakie
CairoMakie.activate!(visible = false)


@time @testset "Unit tests" begin
    include("unit_test.jl")
end

@time @testset "Simulation" begin
    include("integration_test.jl")
end
