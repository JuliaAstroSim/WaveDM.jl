using Test
using WaveDM

@testset "WaveDM Refactoring Verification" begin
    @testset "Module Loading" begin
        @test :WaveDM in names(Main)
    end

    @testset "Core Functions" begin
        @test WavesDM.setup_grid !== nothing  # This should fail - module name issue
    end
end
