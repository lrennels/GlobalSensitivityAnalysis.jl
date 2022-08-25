using Test
using GlobalSensitivityAnalysis

@testset "Test Utils" begin
    @info("test_utils.jl")
    include("test_utils.jl")
end

@testset "Test Sobol Method" begin

    @testset "Compare to SALib.py" begin
        @info("sobol/test_sobol_sobolsequence.jl")
        include("sobol/test_sobol_sobolsequence.jl")

        @info("sobol/test_sobol_ishigami.jl")
        include("sobol/test_sobol_ishigami.jl")

        @info("sobol/test_sobol_borehole.jl")
        include("sobol/test_sobol_borehole.jl")

        @info("sobol/test_sobol_ci.jl")
        include("sobol/test_sobol_ci.jl")
    end

    @testset "Unit Testing" begin
        @info("sobol/test_sobol_unit.jl")
        include("sobol/test_sobol_unit.jl")
    end
end

@testset "Test Delta Method" begin

    @testset "Compare to SALib.py" begin
        include("delta/test_delta_ishigami.jl")
        include("delta/test_delta_borehole.jl")
    end

    @testset "Unit Testing" begin
        include("delta/test_delta_unit.jl")
    end
end
