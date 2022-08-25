module Test_Delta_Borehole

using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures
using GlobalSensitivityAnalysis
# using VegaLite

import GlobalSensitivityAnalysis: borehole, sample
import StatsBase: ordinalrank

################################################################################
## SETUP
################################################################################

# define the (uncertain) parameters of the problem and their distributions
data = DeltaData(
    params = OrderedDict(
        :rw => Normal(0.10, 0.0161812),
        :r => LogNormal(7.71, 1.0056),
        :Tu => Uniform(63070, 115600),
        :Hu => Uniform(990, 1110),
        :Tl => Uniform(63.1, 116),
        :Hl => Uniform(700, 820),
        :L => Uniform(1120, 1680),
        :Kw => Uniform(9855, 12045)
    ),
    N = 1000
)

N = data.N
D = length(data.params)

##
## Sampling 
##

py_samples = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_samples.csv"), header_exists=false, colnames = ["rw", "r", "Tu", "Hu", "Tl", "Hl", "L", "Kw"]) |> DataFrame
julia_samples = sample(data) 

quants = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
sig_level = 0.05
# output_dir = joinpath(@__DIR__, "../output/LHS_QuantileTesting/borehole")
# mkpath(output_dir)

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

##
## Analysis
##

# Get the samples
samples = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_samples.csv"), header_exists=false, colnames = ["rw", "r", "Tu", "Hu", "Tl", "Hl", "L", "Kw"]) |> DataFrame
samples = Matrix(samples)

# check borehole
python_Y = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_borehole.csv"), header_exists=false) |> DataFrame
julia_Y = borehole(samples)
@test python_Y[:,1] ≈ julia_Y atol = ATOL_sample

# julia results
julia_results = analyze(data, samples, julia_Y; num_resamples = 1_000) 

# python results
py_firstorder = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_firstorder.csv"), header_exists=false) |> DataFrame
py_delta = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_delta.csv"), header_exists=false) |> DataFrame
py_firstorder_conf = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_firstorder_conf.csv"), header_exists=false) |> DataFrame
py_delta_conf = load(joinpath(@__DIR__, "../data/delta/py_borehole/py_delta_conf.csv"), header_exists=false) |> DataFrame

# test indices
@test julia_results[:firstorder] ≈ Matrix(py_firstorder) atol = ATOL_delta
@test ordinalrank(julia_results[:firstorder]) == ordinalrank(py_firstorder[!,:Column1])

@test julia_results[:delta] ≈ Matrix(py_delta) atol = 0.05 # TODO - this seems too high?
@test ordinalrank(julia_results[:delta]) == ordinalrank(py_delta[!,:Column1])

# test confidence intervals
@test julia_results[:firstorder_conf] ≈ Matrix(py_firstorder_conf) atol = ATOL_CI
@test julia_results[:delta_conf] ≈ Matrix(py_delta_conf) atol = ATOL_CI

end