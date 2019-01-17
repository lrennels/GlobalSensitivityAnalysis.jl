"""
    SobolParams(names::AbstractArray{String, 1}, dists::AbstractArray{Distribution, 1})

Create a `SobolParams` type containing the `names` of the uncertain parameters
of the problem, and the distributions `dists`, corresponding to each parameter.
"""
mutable struct SobolParams
    names::AbstractArray{String, 1}
    dists::AbstractArray{Distribution, 1}
end

"""
    SobolParams(params::SobolParams, firstorder::AbstractArray{<:Number, 1}, totalorder::AbstractArray{<:Number, 1})

Create a `SobolResults` type containing the `names` of the uncertain parameters
of the problem and the corresponding first order and total order sensitvity indicis 
for each of the parameters, held in `firstorder` and `totalorder` respectively.
"""
mutable struct SobolResults
    params::SobolParams
    firstorder::AbstractArray{<:Number, 1}
    totalorder::AbstractArray{<:Number, 1}
end
