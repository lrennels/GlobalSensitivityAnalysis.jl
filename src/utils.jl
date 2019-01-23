"""
    SobolParams(names::AbstractArray{String, N}, dists::AbstractArray{Distribution, N})

Create a `SobolParams` type containing the `names` of the uncertain parameters
of the problem, and the distributions `dists`, corresponding to each parameter, 
and the `num_samples`.
"""
mutable struct SobolParams
    names::AbstractArray{String, N} where N
    dists::AbstractArray{Distribution, N} where N
    num_samples::Int # commonly referred to as N in the documentation
end

"""
    SobolResults(params::SobolParams, firstorder::AbstractArray{<:Number, N}, totalorder::AbstractArray{<:Number, N}, num_samples::Int)

Create a `SobolResults` type containing the `names` of the uncertain parameters
of the problem and the corresponding first order and total order sensitvity indicis 
for each of the parameters, held in `firstorder` and `totalorder` respectively, 
and the `num_samples`..
"""
mutable struct SobolResults
    params::SobolParams
    firstorder::AbstractArray{<:Number, N} where N
    totalorder::AbstractArray{<:Number, N} where N
    num_samples::Int # commonly referred to as N in the documentation
end
