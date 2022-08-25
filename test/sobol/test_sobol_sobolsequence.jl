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

end