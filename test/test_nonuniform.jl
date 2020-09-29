using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures

import GlobalSensitivityAnalysis: ishigami, split_output

################################################################################
## SET CONSTANTS
################################################################################
ATOL = 1e-9
ATOL_CI = 1e-2

################################################################################
## JULIA
################################################################################

# define the (uncertain) parameters of the problem and their distributions
data = SobolData(
    params = OrderedDict(:x1 => Normal(1, 0.2),
        :x2 => Uniform(0.75, 1.25),
        :x3 => LogNormal(0, 0.5)),
    N = 1000
)

N = data.N
D = length(data.params)

# sampling
julia_samples = sample(data) |> DataFrame
julia_ishigami = ishigami(convert(Matrix, julia_samples)) |> DataFrame

# save this DataFrame for Mimi Testing
# names!(julia_saltelli, Symbol.(params.names))
# save("/Users/lisarennels/JuliaProjects/SAJulia/test/data/nonuniform_samples.csv", julia_saltelli)

# analysis
julia_A, julia_B, julia_AB, julia_BA = split_output(convert(Matrix, julia_ishigami), N, D, data.calc_second_order)
julia_results = analyze(data, convert(Matrix, julia_ishigami), num_resamples = 10_000) 

################################################################################ 
## Python
################################################################################

# sampling
py_samples = load("data/py_nonuniform/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_ishigami = load("data/py_nonuniform/py_ishigami.csv", header_exists=false) |> DataFrame

# analysis
py_A = load("data/py_nonuniform/py_A.csv", header_exists=false) |> DataFrame
py_B = load("data/py_nonuniform/py_B.csv", header_exists=false) |> DataFrame
py_AB = load("data/py_nonuniform/py_AB.csv", header_exists=false) |> DataFrame
py_BA = load("data/py_nonuniform/py_BA.csv", header_exists=false) |> DataFrame

py_firstorder = load("data/py_nonuniform/py_firstorder.csv", header_exists=false) |> DataFrame
py_secondorder = load("data/py_nonuniform/py_secondorder.csv", header_exists=false) |> DataFrame
py_totalorder = load("data/py_nonuniform/py_totalorder.csv", header_exists=false) |> DataFrame

py_firstorder_conf = load("data/py_nonuniform/py_firstorder_conf.csv", header_exists=false) |> DataFrame
py_secondorder_conf = load("data/py_nonuniform/py_secondorder_conf.csv", header_exists=false) |> DataFrame
py_totalorder_conf = load("data/py_nonuniform/py_totalorder_conf.csv", header_exists=false) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Non-Uniform Sampling" begin
    @test convert(Matrix, julia_samples) ≈ convert(Matrix, py_samples) atol = ATOL
end

@testset "Non-Uniform Analysis" begin
    
    @test convert(Matrix, julia_ishigami) ≈ convert(Matrix, py_ishigami) atol = ATOL
    @test julia_A ≈ convert(Matrix, py_A) atol = ATOL
    @test julia_B ≈ convert(Matrix, py_B) atol = ATOL
    @test julia_AB ≈ convert(Matrix, py_AB) atol = ATOL
    @test julia_BA ≈ convert(Matrix, py_BA) atol = ATOL

    @test julia_results[:firstorder] ≈ convert(Matrix, py_firstorder) atol = ATOL
    @test julia_results[:totalorder]≈ convert(Matrix, py_totalorder) atol = ATOL

    for i = 1:D
        for j = i+1:D
            @test julia_results[:secondorder][i,j] ≈ convert(Matrix, py_secondorder)[i,j] atol = ATOL
        end
    end

    # test confidence intervals (using rtol here, big values throw things off)
    @test julia_results[:firstorder_conf] ≈ convert(Matrix, py_firstorder_conf) atol = ATOL_CI
    @test julia_results[:totalorder_conf] ≈ convert(Matrix, py_totalorder_conf) atol = ATOL_CI

    for i = 1:D
        for j = i+1:D
            @test julia_results[:secondorder_conf][i,j] ≈ convert(Matrix, py_secondorder_conf)[i,j] atol = ATOL_CI
        end
    end

end
