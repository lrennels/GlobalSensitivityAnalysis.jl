using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures
# using Plots 

# using VegaLite

import GlobalSensitivityAnalysis: ishigami
import StatsBase: ordinalrank

################################################################################
## SETUP
################################################################################

# define the (uncertain) parameters of the problem and their distributions
data = DeltaData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 1000
)

N = data.N
D = length(data.params)

@testset "Ishigami Function - Sampling" begin

    py_samples = load("data/delta/py_ishigami/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
    julia_samples = GlobalSensitivityAnalysis.sample(data) 

    quants = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    sig_level = 0.05
    output_dir = joinpath(@__DIR__, "../output/LHS_QuantileTesting/ishigami")
    mkpath(output_dir)

    for i = 1:D
        sampleA = py_samples[:,i]
        sampleB = julia_samples[:,i]

        # Make dataframe
        # df = DataFrame(A = sampleA, B = sampleB)
        # save(joinpath(output_dir, "samples$(N)_p$i.csv"), df)

        # # Quick plot comparison
        # p1 = df |> @vlplot(:bar, x = {:A, bin = {step=0.1}}, y="count()")
        # save(joinpath(output_dir, "sampleA_$(N)_p$i.png"), p1)
        # p2 = df |> @vlplot(:bar, x = {:B, bin = {step=0.1}}, y="count()", background=:white)
        # save(joinpath(output_dir, "sampleB_$(N)_p$i.png"), p2)

        # Run WRS quantile matching
        results = WRS.pb2gen(sampleA, sampleB, quantiles = quants) |> DataFrame
        # save(joinpath(output_dir, "LHS Sampling Quantile Comparison N$(N)_p$i.csv"), results)

        @test sum(results[!,:signif]) == 0

    end

end

@testset "Ishigami Function - Analysis" begin

    # Get the samples - here we will use the Python samples to keep things consistent,
    # could also use Julia samples run through both the SALib and GSA functions
    samples = load("data/delta/py_ishigami/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
    samples = Matrix(samples)

    # check ishigami
    python_Y = load("data/delta/py_ishigami/py_ishigami.csv", header_exists=false) |> DataFrame
    julia_Y = ishigami(samples)
    @test python_Y[:,1] ≈ julia_Y atol = ATOL_sample

    # julia results
    julia_results = analyze(data, samples, julia_Y; num_resamples = 1_000) 

    # python results
    py_firstorder = load("data/delta/py_ishigami/py_firstorder.csv", header_exists=false) |> DataFrame
    py_delta = load("data/delta/py_ishigami/py_delta.csv", header_exists=false) |> DataFrame
    py_firstorder_conf = load("data/delta/py_ishigami/py_firstorder_conf.csv", header_exists=false) |> DataFrame
    py_delta_conf = load("data/delta/py_ishigami/py_delta_conf.csv", header_exists=false) |> DataFrame

    # test indices - check values and orderings
    @test julia_results[:firstorder] ≈ Matrix(py_firstorder) atol = ATOL_delta
    @test ordinalrank(julia_results[:firstorder]) == ordinalrank(py_firstorder[!,:Column1])

    @test julia_results[:delta] ≈ Matrix(py_delta) atol = 0.05 # TODO - this seems too high?
    @test ordinalrank(julia_results[:delta]) == ordinalrank(py_delta[!, :Column1])

    # test confidence intervals
    @test julia_results[:firstorder_conf] ≈ Matrix(py_firstorder_conf) atol = ATOL_CI
    @test julia_results[:delta_conf] ≈ Matrix(py_delta_conf) atol = ATOL_CI

end
