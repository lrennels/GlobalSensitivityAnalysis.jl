using Distributions
using DataStructures
using Statistics

"""
    SobolData

A struct which holds all information needed for the sampling and analysis of a
specific problem using Sobol Analysis:

`params::OrderedDict{Symbol, <:Distribution} = Dict(): a dictionary mapping parameter names to their Distribution
`calc_second_order::Bool = false`: whether or not to calculate second order sensitivity indicies
`N::Int = 1000`: the number of runs
`results::Dict{} = Dict()`: the results of the sobol analysis
"""
mutable struct SobolData
    params::OrderedDict{Symbol, <:Distribution} # TODO: error here with non uniform ...
    calc_second_order::Bool
    N::Int 
    results::Dict{}

    function SobolData(params::OrderedDict{Symbol, <:Distribution} = Dict(), calc_second_order::Bool = false, N::Int = 1000, results = Dict())
        return new(params, calc_second_order, N, results)
    end
end


"""
    scale_sobol_seq(sequence::AbstractArray{<:Number, N}, dists::AbstractArray{Distribution, N})

Rescale a Sobol `sequence` of parameters from the 0-to-1 range to their corresponding 
univeariate distributions `dists`.  
"""
function scale_sobol_seq(sequence::AbstractArray{<:Number, N1}, dists::AbstractArray{<:Distribution, N2}) where N1 where N2
    D = length(dists) # number of parameters
    for param in 1:D
        dist = dists[param]
        if !(typeof(dist) <: UnivariateDistribution)
            error("Distribution must be a subtype of UnivariateDistribution")
        end
        sequence[:, [param, param + D]] = quantile.(dist, sequence[:, [param, param + D]])
    end
    return sequence
end
