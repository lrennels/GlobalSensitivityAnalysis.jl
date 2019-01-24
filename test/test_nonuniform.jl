################################################################################
## JULIA
################################################################################

# define the (uncertain) parameters of the problem and their distributions
params = SobolParams(
    ["name1", "name2", "name3"], 
    [Normal(1, 0.2), Uniform(0.75, 1.25), LogNormal(20, 4)],
    100
)

# sampling
julia_sobol = sobol_sequence(params.num_samples, length(params.names)) |> DataFrame 
julia_saltelli = saltelli_sample(params) |> DataFrame
julia_ishigami = ishigami(convert(Matrix, julia_saltelli)) |> DataFrame

# save this DataFrame for Mimi Testing
# names!(julia_saltelli, Symbol.(params.names))
# save("/Users/lisarennels/JuliaProjects/SAJulia/test/data/nonuniform_samples.csv", julia_saltelli)

# analysis
julia_A, julia_B, julia_AB = split_output(convert(Matrix, julia_ishigami), params.num_samples, length(params.names)) 
julia_results = sobol_analyze(params, convert( Matrix, julia_ishigami)) 

################################################################################
## Python
################################################################################

# sampling
py_sobol = load("data/py_nonuniform/py_sobol.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
py_saltelli = load("data/py_nonuniform/py_saltelli.csv", header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame
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
    @test convert(Matrix, julia_sobol) ≈ convert(Matrix, py_sobol) atol = 1e-9
    @test convert(Matrix, julia_saltelli)[:, 2] ≈ convert(Matrix, py_saltelli)[:, 2] atol = 1e-9 # the uniform param
end

@testset "Non-Uniform Analysis" begin
    @test julia_results.firstorder ≈ convert(Matrix, py_firstorder) atol = 1e-9
    @test julia_results.totalorder ≈ convert(Matrix, py_totalorder) atol = 1e-9

end