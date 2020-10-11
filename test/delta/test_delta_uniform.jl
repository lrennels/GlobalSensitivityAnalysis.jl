using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures

import GlobalSensitivityAnalysis: ishigami

################################################################################
## SETUP
################################################################################

ATOL = 1e-9
ATOL_CI = 1e-2

# define the (uncertain) parameters of the problem and their distributions
data = DeltaData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 1000
)

N = data.N
D = length(data.params)

@testset "Uniform Sampling" begin
    # TODO
end


@testset "Uniform Analysis" begin

    # Get the samples
    samples = load("data/delta/py_uniform/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
    
    # check ishigami
    python_Y = load("data/delta/py_uniform/py_ishigami.csv", header_exists=false) |> DataFrame
    julia_Y = ishigami(convert(Array, samples))
    @test python_Y[:1] ≈ julia_Y atol = ATOL

    # julia results
    julia_results = analyze(data, convert(Array, samples), convert(Matrix, Y); num_resamples = 1_000) 

    # python results
    py_firstorder = load("data/delta/py_uniform/py_firstorder.csv", header_exists=false) |> DataFrame
    py_delta = load("data/delta/py_uniform/py_delta.csv", header_exists=false) |> DataFrame
    py_firstorder_conf = load("data/delta/py_uniform/py_firstorder_conf.csv", header_exists=false) |> DataFrame
    py_delta_conf = load("data/delta/py_uniform/py_delta_conf.csv", header_exists=false) |> DataFrame

    # test indices
    @test julia_results[:firstorder] ≈ convert(Matrix, py_firstorder) atol = ATOL
    @test julia_results[:delta] ≈ convert(Matrix, py_delta) atol = ATOL

    # test confidence intervals
    @test julia_results[:firstorder_conf] ≈ convert(Matrix, py_firstorder_conf) atol = ATOL_CI
    @test julia_results[:delta_conf] ≈ convert(Matrix, py_delta_conf) atol = ATOL_CI

end
