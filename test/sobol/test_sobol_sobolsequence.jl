module Test_Sobol_SobolSequence

using Test
using DataFrames
using CSVFiles
using DataStructures
using Sobol
using Distributions
using GlobalSensitivityAnalysis

ATOL_sobol = 1e-9

# define the (uncertain) parameters of the problem and their distributions
# (this data was used for the pre-computed py_sobolseq)
data = SobolData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 1000
)

N = data.N
D = length(data.params)

s = Sobol.SobolSeq(D) 
Sobol_sobolseq = vcat(zeros(1,D), hcat([Sobol.next!(s) for i = 1:N-1]...)') # Sobol package sobol seq function

py_sobolseq = load(joinpath(@__DIR__, "../data/sobol/py_ishigami/py_sobolseq.csv"), header_exists=false, colnames = ["x1", "x2", "x3"]) |> DataFrame

@test convert(Matrix, Sobol_sobolseq)  ≈ Matrix(py_sobolseq)  atol = ATOL_sobol

s = Sobol.SobolSeq(3)
expected = [[0.5, 0.5, 0.5], # note does not start with 0s as do some Sobol sequences
                [0.75, 0.25, 0.25],
                [0.25, 0.75, 0.75],
                [0.375, 0.375, 0.625],
                [0.875, 0.875, 0.125],
                [0.625, 0.125, 0.875],
                [0.125, 0.625, 0.375],
                [0.1875, 0.3125, 0.9375],
                [0.6875, 0.8125, 0.4375]]
for i in 1:9
    @test [Sobol.next!(s)...] == expected[i] # skips teh zeros
end

end