using Test

@testset "Compare to SALib.py" begin
    include("test_sobolsequence.jl")
    include("test_uniform.jl")
    include("test_nonuniform.jl")
end

@testset "Unit Testing" begin
    include("unit_tests.jl")
end
