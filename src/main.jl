using Distributions

include("sobol_sequence.jl")
include("saltelli.jl")
include("test_functions/ishigami.jl")

# define the (uncertain) parameters of the problem and their distributions
params = SobolParams(
    ["x1", "x2", "x3"], 
    [Uniform(-3.14159265359, 3.14159265359), Uniform(-3.14159265359, 3.14159265359), 
    Uniform(-3.14159265359, 3.14159265359)]
)
N = 100000
D = 3

# generate samples using Sobol sequence
param_values = saltelli_sample(params, N)

# generate samples using a given sequence
input_data = rand(2N, D)
rand_param_values = saltelli_sample(params, N, rand_base_sequence = input_data)

# Run model (example)
Y = ishigami(param_values)
