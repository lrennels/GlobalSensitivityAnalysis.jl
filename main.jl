using Distributions

include("sobol_sequence.jl")
include("saltelli.jl")

# define the (uncertain) parameters of the problem and their distributions
params = SobolParams(
    ["x1", "x2", "x3"], 
    [Uniform(-3.14159265359, 3.14159265359), Uniform(-3.14159265359, 3.14159265359), 
    Uniform(-3.14159265359, 3.14159265359)]
)

# generate samples
saltelli_sample(params, 5)
