using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures

import GlobalSensitivityAnalysis: borehole, split_output, sample

################################################################################
## JULIA
################################################################################

# define the (uncertain) parameters of the problem and their distributions
data = SobolData(
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

# sampling
julia_samples = sample(data) |> DataFrame
julia_borehole = borehole(convert(Matrix, julia_samples)) |> DataFrame

# analysis
julia_A, julia_B, julia_AB, julia_BA = split_output(convert(Matrix, julia_borehole), N, D, data.calc_second_order)
julia_results = analyze(data, convert(Matrix, julia_borehole), num_resamples = 10_000) 

################################################################################ 
## Python
################################################################################

# sampling
py_samples = load("data/sobol/py_borehole/py_samples.csv", header_exists=false, colnames = ["rw", "r", "Tu", "Hu", "Tl", "Hl", "L", "Kw"]) |> DataFrame
py_borehole = load("data/sobol/py_borehole/py_borehole.csv", header_exists=false) |> DataFrame

# analysis
py_A = load("data/sobol/py_borehole/py_A.csv", header_exists=false) |> DataFrame
py_B = load("data/sobol/py_borehole/py_B.csv", header_exists=false) |> DataFrame
py_AB = load("data/sobol/py_borehole/py_AB.csv", header_exists=false) |> DataFrame
py_BA = load("data/sobol/py_borehole/py_BA.csv", header_exists=false) |> DataFrame

py_firstorder = load("data/sobol/py_borehole/py_firstorder.csv", header_exists=false) |> DataFrame
py_secondorder = load("data/sobol/py_borehole/py_secondorder.csv", header_exists=false) |> DataFrame
py_totalorder = load("data/sobol/py_borehole/py_totalorder.csv", header_exists=false) |> DataFrame

py_firstorder_conf = load("data/sobol/py_borehole/py_firstorder_conf.csv", header_exists=false) |> DataFrame
py_secondorder_conf = load("data/sobol/py_borehole/py_secondorder_conf.csv", header_exists=false) |> DataFrame
py_totalorder_conf = load("data/sobol/py_borehole/py_totalorder_conf.csv", header_exists=false) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Borehole Function - Sampling" begin
    @test convert(Matrix, julia_samples) ≈ convert(Matrix, py_samples) atol = ATOL_sample
end

@testset "Borehole Function - Analysis" begin
    
    @test convert(Matrix, julia_borehole) ≈ convert(Matrix, py_borehole) atol = ATOL_sobol
    @test julia_A ≈ convert(Matrix, py_A) atol = ATOL_sobol
    @test julia_B ≈ convert(Matrix, py_B) atol = ATOL_sobol
    @test julia_AB ≈ convert(Matrix, py_AB) atol = ATOL_sobol
    @test julia_BA ≈ convert(Matrix, py_BA) atol = ATOL_sobol

    @test julia_results[:firstorder] ≈ convert(Matrix, py_firstorder) atol = ATOL_sobol
    @test julia_results[:totalorder]≈ convert(Matrix, py_totalorder) atol = ATOL_sobol

    for i = 1:D
        for j = i+1:D
            @test julia_results[:secondorder][i,j] ≈ convert(Matrix, py_secondorder)[i,j] atol = ATOL_sobol
        end
    end

    # test confidence intervals
    @test julia_results[:firstorder_conf] ≈ convert(Matrix, py_firstorder_conf) atol = ATOL_CI
    @test julia_results[:totalorder_conf] ≈ convert(Matrix, py_totalorder_conf) atol = ATOL_CI

    for i = 1:D
        for j = i+1:D
            @test julia_results[:secondorder_conf][i,j] ≈ convert(Matrix, py_secondorder_conf)[i,j] atol = ATOL_CI
        end
    end

end
