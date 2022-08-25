module Test_Utils

using Test
using Distributions
using GlobalSensitivityAnalysis

# scale_samples!
seq = rand(100, 8)
original_seq = copy(seq)
dists = [Normal(1, 0.2), Uniform(0.75, 1.25), LogNormal(0, 0.5), TriangularDist(0, 4, 1)]
GlobalSensitivityAnalysis.scale_samples!(seq, dists)
@test size(seq) == size(original_seq)
@test seq != original_seq

end