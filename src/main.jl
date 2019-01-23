using Distributions

include("sobol_sequence.jl")
include("saltelli.jl")
include("sobol_analysis.jl")
include("test_functions/ishigami.jl")

# define the (uncertain) parameters of the problem and their distributions as well
# as the number of samples
params = SobolParams(
    ["x1", "x2", "x3"], 
    [Uniform(-3.14159265359, 3.14159265359), Uniform(-3.14159265359, 3.14159265359), 
    Uniform(-3.14159265359, 3.14159265359)],
    1000
)

# generate samples using Sobol sequence
param_values = saltelli_sample(params)

# run model (example)
Y = ishigami(param_values)

# perform Sobol Analysis
results = sobol_analyze(params, Y)
