using Test
using GlobalSensitivityAnalysis

@testset "Test Utils" begin
    include("test_utils.jl")
end

@testset "Test Sobol Method" begin

    ATOL = 1e-9
    ATOL_CI = 1e-2

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

    ATOL = 1e-9
    ATOL_IDX = 1e-1
    ATOL_CI = 1e-2
    ATOL_SAMPLE = 1e-2

    @testset "Compare to SALib.py" begin
        include("delta/test_delta_uniform.jl")
        include("delta/test_delta_nonuniform.jl")
        include("delta/test_delta_ci.jl")
    end

    @testset "Unit Testing" begin
        include("delta/test_delta_unit.jl")
    end
end
