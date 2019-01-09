using Distributions
using Test
using DataFrames
using CSVFiles

################################################################################
## JULIA
################################################################################
include("../src/sobol_sequence.jl")
include("../src/saltelli.jl")
include("../src/test_functions/ishigami.jl")

# define number of samples
N = 100000
D = 3

# define the (uncertain) parameters of the problem and their distributions
params = SobolParams(
    ["x1", "x2", "x3"], 
    [Uniform(-3.14159265359, 3.14159265359), Uniform(-3.14159265359, 3.14159265359), 
    Uniform(-3.14159265359, 3.14159265359)]
)

# get arrays for comparison
julia_sobol = sobol_sequence(N, D) |> DataFrame
julia_saltelli = saltelli_sample(params, N) |> DataFrame
julia_ishigami = ishigami(convert(Matrix, julia_saltelli)) |> DataFrame

################################################################################
## Python
################################################################################

py_sobol = load("data/py_sobol.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_saltelli = load("data/py_saltelli.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_ishigami = load("data/py_ishigami.csv", header_exists=false) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Sampling" begin
    @test convert(Matrix, julia_sobol) ≈ convert(Matrix, py_sobol) atol = 1e-9
    @test convert(Matrix, julia_saltelli) ≈ convert(Matrix, py_saltelli) atol = 1e-9
end

@testset "Analysis" begin
    @test convert(Matrix, julia_ishigami) ≈ convert(Matrix, py_ishigami) atol = 1e-9
end
