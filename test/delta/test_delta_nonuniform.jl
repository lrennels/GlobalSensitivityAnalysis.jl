using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures

import GlobalSensitivityAnalysis: ishigami

################################################################################
## SETUP
################################################################################

# define the (uncertain) parameters of the problem and their distributions
data = DeltaData(
    params = OrderedDict(:x1 => Normal(1, 0.2),
        :x2 => Uniform(0.75, 1.25),
        :x3 => LogNormal(0, 0.5)),
    N = 1000
)

N = data.N
D = length(data.params)

# TODO: get the python files from Jupyter

# TODO is this the right way to compare the samples?
# see: https://www.nature.com/articles/sdata2018187#Sec5

@testset "Non-Uniform Sampling" begin

    py_samples = convert(Matrix, load("data/delta/py_nonuniform/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame)
    julia_samples = sample(data) 

    for i = 1:D
        julia_quantiles = quantile(julia_samples[:,i])
        py_quantiles = quantile(py_samples[:,i])
        @test julia_quantiles ≈ py_quantiles atol = ATOL_SAMPLE
    end

end


@testset "Non-Uniform Analysis" begin

    # Get the samples
    samples = load("data/delta/py_nonuniform/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
    
    # check ishigami
    python_Y = load("data/delta/py_nonuniform/py_ishigami.csv", header_exists=false) |> DataFrame
    julia_Y = ishigami(convert(Array, samples))
    @test python_Y[:,1] ≈ julia_Y atol = ATOL

    # julia results
    julia_results = analyze(data, convert(Array, samples), julia_Y; num_resamples = 1_000) 

    # python results
    py_firstorder = load("data/delta/py_nonuniform/py_firstorder.csv", header_exists=false) |> DataFrame
    py_delta = load("data/delta/py_nonuniform/py_delta.csv", header_exists=false) |> DataFrame
    py_firstorder_conf = load("data/delta/py_nonuniform/py_firstorder_conf.csv", header_exists=false) |> DataFrame
    py_delta_conf = load("data/delta/py_nonuniform/py_delta_conf.csv", header_exists=false) |> DataFrame

    # test indices
    @test julia_results[:firstorder] ≈ convert(Matrix, py_firstorder) atol = ATOL_IDX
    @test julia_results[:delta] ≈ convert(Matrix, py_delta) atol = ATOL_IDX

    # test confidence intervals
    @test julia_results[:firstorder_conf] ≈ convert(Matrix, py_firstorder_conf) atol = ATOL_CI
    @test julia_results[:delta_conf] ≈ convert(Matrix, py_delta_conf) atol = ATOL_CI

end