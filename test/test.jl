using Distributions
using Test
using DataFrames

################################################################################
## JULIA
################################################################################
include("../src/sobol_sequence.jl")
include("../src/saltelli.jl")

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

################################################################################
## Python
################################################################################

py_sobol = load("data/py_sobol.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_saltelli = load("data/py_saltelli.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Sampling" begin
    @test convert(Matrix, julia_sobol) ≈ convert(Matrix, py_sobol) atol = 1e-9

    # check pieces of py_saltelli
    julia_A = julia_saltelli[1:N, :]
    py_A = py_saltelli[1:N, :]
    @test convert(Matrix, julia_A) ≈ convert(Matrix, py_A) atol = 1e-3

    # check pieces of py_saltelli (A, ABs, and B) because currently the Python
    # version from SALib also includes the BAs for second order, but we do not 
    # include those
    julia_A = julia_saltelli[1:N, :]
    py_A = py_saltelli[1:N, :]
    @test convert(Matrix, julia_A) ≈ convert(Matrix, py_A) atol = 1e-3

    julia_ABs = julia_saltelli[N+1:end - N, :]
    py_ABs = py_saltelli[N+1:(N+1) + (N*D) - 1, :]
    @test convert(Matrix, julia_A) ≈ convert(Matrix, py_A) atol = 1e-3

    julia_B = julia_saltelli[end-N+1:end, :]
    py_B = py_saltelli[end-N+1:end, :]
    @test convert(Matrix, julia_B) ≈ convert(Matrix, py_B) atol = 1e-3

end