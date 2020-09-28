using Test

@testset "Test Utils" begin
    include("test_utils.jl")
end

@testset "Test Sobol Method" begin

    @testset "Compare to SALib.py" begin
        include("sobol/test_sobol_sobolsequence.jl")
        include("sobol/test_sobol_uniform.jl")
        include("sobol/test_sobol_nonuniform.jl")
        include("sobol/test_sobol_ci.jl")
    end

    @testset "Unit Testing" begin
        include("sobol/test_sobol_unit.jl")
    end
end

@testset "Test Delta Method" begin

    @testset "Compare to SALib.py" begin
        # TODO
    end

    @testset "Unit Testing" begin
        include("delta/test_delta_unit.jl")
    end
end
