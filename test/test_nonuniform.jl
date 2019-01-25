using Distributions
using Test
using DataFrames
using CSVFiles
using DataStructures

################################################################################
## JULIA
################################################################################
include("../src/sample_sobol.jl")
include("../src/analyze_sobol.jl")
include("../src/test_functions/ishigami.jl")

# define the (uncertain) parameters of the problem and their distributions
data = SobolData(
    params = OrderedDict(:x1 => Normal(1, 0.2),
        :x2 => Uniform(0.75, 1.25),
        :x3 => LogNormal(20, 4)),
    N = 100
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
julia_A, julia_B, julia_AB = split_output(convert(Matrix, julia_ishigami), N, D)
julia_results = analyze(data, convert( Matrix, julia_ishigami)) 

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
py_firstorder = load("data/py_nonuniform/py_firstorder.csv", header_exists=false) |> DataFrame
py_totalorder = load("data/py_nonuniform/py_totalorder.csv", header_exists=false) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Non-Uniform Sampling" begin
    @test convert(Matrix, julia_samples)[:, 1:2] ≈ convert(Matrix, py_samples)[:, 1:2] atol = 1e-9 # lognormal is different
end

@testset "Non-Uniform Analysis" begin
    @test julia_results[:firstorder] ≈ convert(Matrix, py_firstorder) atol = 1e-9
    @test julia_results[:totalorder]≈ convert(Matrix, py_totalorder) atol = 1e-9
end
