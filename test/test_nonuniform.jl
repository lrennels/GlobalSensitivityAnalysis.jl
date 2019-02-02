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
        :x3 => LogNormal(20, 4),
        :x4 => TriangularDist(0, 4, 1)),
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
julia_A, julia_B, julia_AB, julia_BA = split_output(convert(Matrix, julia_ishigami), N, D, data.calc_second_order)
julia_results = analyze(data, convert( Matrix, julia_ishigami)) 

################################################################################
## Python
################################################################################

# sampling
py_samples = load("data/py_nonuniform/py_samples.csv", header_exists=false, colnames = ["x1", "x2", "x3", "x4"]) |> DataFrame
py_ishigami = load("data/py_nonuniform/py_ishigami.csv", header_exists=false) |> DataFrame

# analysis
py_A = load("data/py_nonuniform/py_A.csv", header_exists=false) |> DataFrame
py_B = load("data/py_nonuniform/py_B.csv", header_exists=false) |> DataFrame
py_AB = load("data/py_nonuniform/py_AB.csv", header_exists=false) |> DataFrame
py_BA = load("data/py_nonuniform/py_BA.csv", header_exists=false) |> DataFrame
py_firstorder = load("data/py_nonuniform/py_firstorder.csv", header_exists=false) |> DataFrame
py_secondorder = load("data/py_nonuniform/py_secondorder.csv", header_exists=false) |> DataFrame
py_totalorder = load("data/py_nonuniform/py_totalorder.csv", header_exists=false) |> DataFrame

################################################################################
## Testing
################################################################################

@testset "Non-Uniform Sampling" begin
    for i = 1:D
        if i != 3 #lognormal treated differently because big values
            @test convert(Matrix, julia_samples)[:, i] ≈ convert(Matrix, py_samples)[:, i] atol = 1e-9 # lognormal numbers too big for atol
        else
            # lognormal scales to huge numbers, so need to test differently
            julia_lognorm = convert(Matrix, julia_samples)[:,i]
            py_lognorm = convert(Matrix, py_samples)[:,i]
            @test maximum((julia_lognorm .- py_lognorm) ./ py_lognorm) < 1e-9
        end
    end
end

@testset "Non-Uniform Analysis" begin
    @test julia_results[:firstorder] ≈ convert(Matrix, py_firstorder) atol = 1e-9
    @test julia_results[:totalorder]≈ convert(Matrix, py_totalorder) atol = 1e-9

    @test julia_results[:totalorder] ≈ convert(Matrix, py_totalorder) atol = 1e-9

    for i = 1:D
        for j = i+1:D
            @test julia_results[:secondorder][i,j] ≈ convert(Matrix, py_secondorder)[i,j] atol = 1e-9
        end
    end
end
