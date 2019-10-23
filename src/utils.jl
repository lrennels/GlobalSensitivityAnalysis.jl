using Distributions
using DataStructures
using Statistics

"""
    SobolData

A struct which holds all information needed for the sampling and analysis of a
specific problem using Sobol Analysis:

`params::Union{OrderedDict{Symbol, <:Any}, Nothing} = nothing`: a dictionary mapping parameter names to their Distribution
`calc_second_order::Bool = false`: whether or not to calculate second order sensitivity indices
`N::Int = 1000`: the number of runs
"""
mutable struct SobolData
    params::Union{OrderedDict{Symbol, <:Any}, Nothing}
    calc_second_order::Bool
    N::Int 

    function SobolData(;params= nothing, calc_second_order = true, N = 1000)
        return new(params, calc_second_order, N)
    end
end


"""
    scale_sobol_seq!(sequence::AbstractArray{<:Number, N}, dists::AbstractArray{T, N})

Rescale a Sobol `sequence` of parameters from the 0-to-1 range to their corresponding 
univeariate distributions `dists`.  
"""
function scale_sobol_seq!(sequence::AbstractArray{<:Number, N1}, dists::AbstractArray{T, N2}) where T where N1 where N2
    D = length(dists) # number of parameters
    for param in 1:D
        dist = dists[param]
        if typeof(dist) <: UnivariateDistribution
            sequence[:, [param, param + D]] = quantile.(dist, sequence[:, [param, param + D]])
        # temporary fix to problem calling quantile. on a EmpiricalDistribution from Mimi
        else
            sequence[:, [param, param + D]] = quantile(dist, sequence[:, [param, param + D]])
        end
    end
end
