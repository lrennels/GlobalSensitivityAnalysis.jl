module GlobalSensitivityAnalysis

using DataStructures
using Distributions
using Statistics
using Sobol

export 
    sample, 
    analyze, 
    SobolData, 
    DeltaData

include("utils.jl")
include("sample_sobol.jl")
include("sample_delta.jl")
include("analyze_sobol.jl")
include("analyze_delta.jl")
include("test_functions/ishigami.jl")
include("test_functions/borehole.jl")

end # module
