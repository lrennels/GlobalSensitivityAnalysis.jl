using Test

@testset "SALib" begin
    include("test_sobolsequence.jl")
    include("test_uniform.jl")
    include("test_nonuniform.jl")
end
