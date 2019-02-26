module GlobalSensitivityAnalysis

using DataStructures
using Distributions
using Statistics
using Sobol

export
    sample, analyze, SobolData

include("utils.jl")
include("sample_sobol.jl")
include("analyze_sobol.jl")
include("test_functions/ishigami.jl")

end # module
