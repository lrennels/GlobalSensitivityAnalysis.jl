#=
Code adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
library for sensitivity analysis. Journal of Open Source Software, 2(9)
=#

"""
    ishigami(param_values::AbstractArray{<:Number, N})

Return an `N` by 1 array of results of the non-monotonic Ishigami Function using the
`N` by 3 array of `param_values`.
"""
function ishigami(param_values::AbstractArray{<:Number, N}) where N

    numruns = size(param_values, 1)

    Y = zeros(numruns, 1)
    A = 7
    B = 0.1

    for i in 1:numruns
        Y[i] = sin(param_values[i, 1]) + A * (sin(param_values[i, 2])^2) + B * (param_values[i, 3] ^ 4) * sin(param_values[i, 1])
    end

    return Y
end
