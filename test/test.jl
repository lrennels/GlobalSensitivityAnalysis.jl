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
include("../src/sobol_analysis.jl")

# define number of samples
N = 100000
D = 3

# define the (uncertain) parameters of the problem and their distributions
params = SobolParams(
    ["x1", "x2", "x3"], 
    [Uniform(-3.14159265359, 3.14159265359), Uniform(-3.14159265359, 3.14159265359), 
    Uniform(-3.14159265359, 3.14159265359)]
)

# sampling
julia_sobol = sobol_sequence(N, D) |> DataFrame
julia_saltelli = saltelli_sample(params, N) |> DataFrame
julia_ishigami = ishigami(convert(Matrix, julia_saltelli)) |> DataFrame

# analysis
julia_A, julia_B, julia_AB = split_output(convert(Matrix, julia_ishigami), N, D)
julia_sobol = sobol_analyze(convert(Matrix, julia_ishigami))

################################################################################
## Python
################################################################################

# sampling
py_sobol = load("data/py_sobol.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_saltelli = load("data/py_saltelli.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_ishigami = load("data/py_ishigami.csv", header_exists=false) |> DataFrame

# analysis
py_A = load("data/py_A.csv", header_exists=false) |> DataFrame
py_B = load("data/py_B.csv", header_exists=false) |> DataFrame
py_AB = load("data/py_AB.csv", header_exists=false) |> DataFrame
py_firstorder = load("data/py_firstorder.csv", header_exists=false) |> DataFrame
py_totalorder = load("data/py_totalorder.csv", header_exists=false) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Sampling" begin
    @test convert(Matrix, julia_sobol) ≈ convert(Matrix, py_sobol) atol = 1e-9
    @test convert(Matrix, julia_saltelli) ≈ convert(Matrix, py_saltelli) atol = 1e-9
end

@testset "Analysis" begin
    @test convert(Matrix, julia_ishigami) ≈ convert(Matrix, py_ishigami) atol = 1e-9
    @test convert(Matrix, julia_A) ≈ convert(Matrix, py_A) atol = 1e-9
    @test convert(Matrix, julia_B) ≈ convert(Matrix, py_B) atol = 1e-9
    @test convert(Matrix, julia_AB) ≈ convert(Matrix, py_AB) atol = 1e-9
    @test julia_sobol.firstorder ≈ py_firstorder atol = 1e-9
    @test julia_sobol.totalorder ≈ py_totalorder atol = 1e-9

end
