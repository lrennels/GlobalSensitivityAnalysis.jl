module GlobalSensitivityAnalysis

using DataStructures
using Distributions
using Statistics
using Sobol

export
    sample, analyze

include("utils.jl")
include("sample_sobol.jl")
include("analyze_sobol.jl")
# include("sobol_sequence.jl") # removed because we are using Sobol.jl for now
include("test_functions/ishigami.jl")

end # module
