module Test_Delta_Regression

using Distributions
using Test
using DataStructures
using GlobalSensitivityAnalysis

import GlobalSensitivityAnalysis: ishigami, sample

# define the (uncertain) parameters of the problem and their distributions
data = DeltaData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 10_000
)

N = data.N
D = length(data.params)

samples = sample(data) 
Y = ishigami(samples)

# test for regression - tolerances as indicated by SALib
results = analyze(data, samples, Y; num_resamples = 10)
@test isapprox(results[:delta], [0.210, 0.358, 0.155], atol=1e-2, rtol=1e-1)
@test isapprox(results[:firstorder], [0.31, 0.44, 0.00], atol=1e-2, rtol=1e-1)

end
