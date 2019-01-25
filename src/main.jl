using Distributions
using DataStructures

include("sample_sobol.jl")
include("analyze_sobol.jl")
include("test_functions/ishigami.jl")

# define the data
data = SobolData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 1000
)

# generate samples using Sobol sequence
samples = sample(data)

# run model (example)
Y = ishigami(samples)

# perform Sobol Analysis
analyze(data, Y)
