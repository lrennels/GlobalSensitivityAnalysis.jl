using Distributions
using DataStructures
using Statistics

"""
    SobolData

A struct which holds all information needed for the sampling and analysis of a
specific problem using Sobol Analysis:

`params::Union{OrderedDict{Symbol, <:Any}, Nothing} = nothing`: a dictionary mapping parameter names to their Distribution
`calc_second_order::Bool = true`: whether or not to calculate second order sensitivity indices
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
    DeltaData

A struct which holds all information needed for the sampling and analysis of a
specific problem using Delta Analysis:

`params::Union{OrderedDict{Symbol, <:Any}, Nothing} = nothing`: a dictionary mapping parameter names to their Distribution
`N::Int = 1000`: the number of runs
"""
mutable struct DeltaData
    params::Union{OrderedDict{Symbol, <:Any}, Nothing}
    N::Int 

    function DeltaData(;params= nothing, N = 1000)
        return new(params, N)
    end
end

"""
    scale_samples!(sequence::AbstractArray{<:Number, N1}, dists::AbstractArray{T, N2})

Rescale a Sobol `sequence` of parameters from the 0-to-1 range to their corresponding 
univariate distributions `dists`. 
"""
function scale_samples!(sequence::AbstractArray{<:Number, N1}, dists::AbstractArray{T, N2}) where T where N1 where N2
    D = length(dists) # number of parameters
    duplicate = size(sequence,2) == 2*D # sobol requires twice as many columns as D

    for param in 1:D
        dist = dists[param]
        if duplicate
            if typeof(dist) <: UnivariateDistribution
                sequence[:, [param, param + D]] = quantile.(dist, sequence[:, [param, param + D]])
            # temporary fix to problem calling quantile on a EmpiricalDistribution from Mimi
            else
                sequence[:, [param, param + D]] = quantile(dist, sequence[:, [param, param + D]])
            end
        else
            if typeof(dist) <: UnivariateDistribution
                sequence[:, param] = quantile.(dist, sequence[:, param])
            # temporary fix to problem calling quantile on a EmpiricalDistribution from Mimi
            else
                sequence[:, param] = quantile(dist, sequence[:, param])
            end
        end
    end
end

"""
    _check_conf_flag(num_resamples::Union{Nothing, Int}, conf_level::Union{Nothing, Number})        
Check to see if confdience interval should be calculated based on the provided
`num_resamples` and `conf_level`.  Error if only one has a value and the other is Nothing.
"""
function _check_conf_flag(num_resamples::Union{Nothing, Int}, conf_level::Union{Nothing, Number})
    num_nothings = (num_resamples === nothing) + (conf_level === nothing)
    if num_nothings == 1
        error("Number of resamples is $num_resamples, while confidence level is $conf_level ... either none or both must be nothing")
    elseif num_nothings == 2
        conf_flag = false
    else
        conf_flag = true
    end
    return conf_flag
end
