using Test
using GlobalSensitivityAnalysis

ATOL_sample = 1e-9
ATOL_sobol = 1e-9 # sobol method indices
ATOL_delta = 1e-3 # delta method indices
ATOL_CI = 1e-2

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
        include("delta/test_delta_uniform.jl")
        include("delta/test_delta_nonuniform.jl")
    end

    @testset "Unit Testing" begin
        include("delta/test_delta_unit.jl")
    end
end
