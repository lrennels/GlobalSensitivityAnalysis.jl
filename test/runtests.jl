using Test

@testset "Test Delta Method" begin
    # TODO
end 

@testset "Test Sobol Method" begin
    @testset "Compare to SALib.py" begin
        include("sobol/test_sobolsequence.jl")
        include("sobol/test_uniform.jl")
        include("sobol/test_nonuniform.jl")
    end

    @testset "Unit Testing" begin
        include("sobol/unit_tests.jl")
    end
end
