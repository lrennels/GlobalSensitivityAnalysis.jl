using Statistics

include("utils.jl")

#=
Code adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
library for sensitivity analysis. Journal of Open Source Software, 2(9)

References
----------
    [1] Sobol, I. M. (2001).  "Global sensitivity indices for nonlinear
        mathematical models and their Monte Carlo estimates."  Mathematics
        and Computers in Simulation, 55(1-3):271-280,
        doi:10.1016/S0378-4754(00)00270-6.
    [2] Saltelli, A. (2002).  "Making best use of model evaluations to
        compute sensitivity indices."  Computer Physics Communications,
        145(2):280-297, doi:10.1016/S0010-4655(02)00280-1.
    [3] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
        S. Tarantola (2010).  "Variance based sensitivity analysis of model
        output.  Design and estimator for the total sensitivity index."
        Computer Physics Communications, 181(2):259-270,
        doi:10.1016/j.cpc.2009.09.018.
=#

"""
    sobol_analyze(params::SobolParams, model_output::AbstractArray{<:Number, 2})

Performs a Sobol Analysis on the `model_output` produced using samples defined 
by the parameters in `params` and returns a SobolResults struct holding the parameters
as well as the first order and total order sensitvity indicies for each of the
uncertain parameters.
"""
# TODO - include second order effects
# TODO - include groups
function sobol_analyze(params::SobolParams, model_output::AbstractArray{<:Number, 2})

    # define constants
    D = length(params.names)
    N = Int(size(model_output, 1) / (D + 2))

    # normalize model output
    model_output = (model_output .- mean(model_output)) ./ std(model_output)

    # separate the model_output into results from matrices "A". "B" and "AB" 
    A, B, AB = split_output(model_output, N, D)

    # compute indicies and produce results
    firstorder = Array{Float64}(undef, D)
    totalorder = Array{Float64}(undef, D)

    for i in 1:D
        firstorder[i] = first_order(A, AB[:, i], B)
        totalorder[i] = total_order(A, AB[:, i], B)
    end

    results = SobolResults(params, firstorder, totalorder)
    
    return results

end

"""
    first_order(A::AbstractArray{<:Number, 1}, AB::AbstractArray{<:Number, 2}, B::AbstractArray{<:Number, 1})

Calculate the first order sensitivity indicies for model outputs given model outputs
separated out into `A`, `AB`, and `A` and normalize by the variance of `[A B]`. [Saltelli et al., 
2010 Table 2 eq (b)]
"""
function first_order(A::AbstractArray{<:Number, 1}, AB::AbstractArray{<:Number, 1}, B::AbstractArray{<:Number, 1})
    return (mean(B .* (AB .- A), dims = 1) / var(vcat(A, B), corrected = false))[1]
end

"""
    total_order(A::AbstractArray{<:Number, 1}, AB::AbstractArray{<:Number, 2}, B::AbstractArray{<:Number, 1})

Calculate the total order sensitivity indicies for model outputs given model outputs
separated out into `A`, `AB`, and `A` and normalize by the variance of `[A B]`. [Saltelli et al., 
2010 Table 2 eq (f)].
"""
function total_order(A::AbstractArray{<:Number, 1}, AB::AbstractArray{<:Number, 1}, B::AbstractArray{<:Number, 1})
    return (0.5 * mean((A .- AB).^2, dims = 1) / var(vcat(A, B), corrected = false))[1]
end

"""
    split_output(model_output::AbstractArray{<:Number, 2}, N, D)

Separate the `model_outputs` into matrices "A", "B", and "AB" for calculation of sensitvity 
indices and return those three matrices.
"""
function split_output(model_output::AbstractArray{<:Number, 2}, N, D)
    stepsize = D + 2

    A = model_output[1:stepsize:end]
    B = model_output[stepsize:stepsize:end]
    
    AB = Array{Float64}(undef, N, D)
    for i in 1:D
        AB[:, i] = model_output[i+1:stepsize:end, :]
    end

    return A, B, AB
end
