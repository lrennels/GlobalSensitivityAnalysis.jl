using Distributions
using PyCall
using Test

################################################################################
## JULIA
################################################################################
include("sobol_sequence.jl")
include("saltelli.jl")

# define number of samples
N = 100000

# define the (uncertain) parameters of the problem and their distributions
params = SobolParams(
    ["x1", "x2", "x3"], 
    [Uniform(-3.14159265359, 3.14159265359), Uniform(-3.14159265359, 3.14159265359), 
    Uniform(-3.14159265359, 3.14159265359)]
)

# get matrices
julia_sobol = sobol_sequence(N, 10)
julia_saltelli = saltelli_sample(params, N)

################################################################################
## Python
################################################################################

Python"""


"""
################################################################################
## Testing
################################################################################

@testset "Sampling" begin

end