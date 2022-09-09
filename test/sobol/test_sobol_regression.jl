module Test_Sobol_Regression

using Test
using GlobalSensitivityAnalysis
using Distributions
using DataStructures

import GlobalSensitivityAnalysis: ishigami, split_output, sample

# define the (uncertain) parameters of the problem and their distributions
data = SobolData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 10_000
)

N = data.N
D = length(data.params)

# sampling
samples = sample(data)
Y = ishigami(samples)

# analysis
results = analyze(data, Y; num_resamples = 10) 

# test for regression - tolerances as indicated by SALib
@test isapprox(results[:firstorder], [0.31, 0.44, 0.00], atol=1e-2, rtol=1e-1)
@test isapprox([results[:secondorder][1,2], results[:secondorder][1,3], results[:secondorder][2,3]], [0.00, 0.25, 0.00], atol=1e-2, rtol=1e-1)
@test isapprox(results[:totalorder], [0.55, 0.44, 0.24], atol=1e-2, rtol=1e-1)

end