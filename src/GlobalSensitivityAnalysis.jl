module GlobalSensitivityAnalysis

using DataStructures
using Distributed
using Distributions
using KernelDensity
using ProgressMeter
using Random
using Sobol
using Statistics 
using StatsBase
using Trapz # alternative considered was NumericalIntegration

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
