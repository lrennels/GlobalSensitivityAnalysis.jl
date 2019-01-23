"""
    SobolParams(names::AbstractArray{String, 1}, dists::AbstractArray{Distribution, 1})

Create a `SobolParams` type containing the `names` of the uncertain parameters
of the problem, and the distributions `dists`, corresponding to each parameter, 
and the `num_samples`.
"""
mutable struct SobolParams
    names::AbstractArray{String, 1}
    dists::AbstractArray{Distribution, 1}
    num_samples::Int # commonly referred to as N in the documentation
end

"""
    SobolParams(params::SobolParams, firstorder::AbstractArray{<:Number, 1}, totalorder::AbstractArray{<:Number, 1})

Create a `SobolResults` type containing the `names` of the uncertain parameters
of the problem and the corresponding first order and total order sensitvity indicis 
for each of the parameters, held in `firstorder` and `totalorder` respectively, 
and the `num_samples`..
"""
mutable struct SobolResults
    params::SobolParams
    firstorder::AbstractArray{<:Number, 1}
    totalorder::AbstractArray{<:Number, 1}
    num_samples::Int # commonly referred to as N in the documentation
end
