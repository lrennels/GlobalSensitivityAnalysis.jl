module SALib

using DataStructures
using Distributions
using Statistics

export
    SobolData, sobol_sequence, sample, analyze!, ishigami, split_output
    
include("utils.jl")

# After loading types, the rest can just be alphabetical
include("analyze_sobol.jl")
include("directions.jl")
include("sample_sobol.jl")
include("sobol_sequence.jl")
include("test_functions/ishigami.jl")

end # module
