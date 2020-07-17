using Test

@testset "Test Sobol Method" begin

    @testset "Compare to SALib.py" begin
        include("sobol/test_sobolsequence.jl")
        include("sobol/test_uniform.jl")
        include("sobol/test_nonuniform.jl")
    end

    @testset "Unit Testing" begin
        include("sobol/unit_tests.jl")
    end

@testset "Test Delta Method" begin

    @testset "Compare to SALib.py" begin
        # TODO
    end

    @testset "Unit Testing" begin
        # TODO
    end

end

end
