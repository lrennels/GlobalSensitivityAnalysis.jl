# Adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
# library for sensitivity analysis. Journal of Open Source Software, 2(9)

"""
    borehole(param_values::AbstractArray{<:Number, N})

Return an `N` by 1 array of results of the Borehole function (http://www.sfu.ca/~ssurjano/borehole.html)
using the `N` by 8 array of `param_values`.
"""
function borehole(param_values::AbstractArray{<:Number, N}) where N

    numruns = size(param_values, 1)
    Y = zeros(numruns, 1)

    for i in 1:numruns

        # parameters
        rw =    param_values[i, 1] # radius of borehole (m)
        r =     param_values[i, 2] # radius of influence (m)
        Tu =    param_values[i, 3] # transmissivity of upper aquifer (m2/yr)
        Hu =    param_values[i, 4] # potentiometric head of upper aquifer (m)
        Tl =    param_values[i, 5] # transmissivity of lower aquifer (m2/yr)
        Hl =    param_values[i, 6] # potentiometric head of lower aquifer (m)
        L =     param_values[i, 7] # length of borehole (m)
        Kw =    param_values[i, 8] # hydraulic conductivity of borehole (m/yr)

        # calculation
        l = log(r/rw)
        num = 2 * pi * Tu * (Hu - Hl)
        denom = l * (1 + ((2 * L * Tu) / (l * rw^2 * Kw)) + (Tu/Tl))
        Y[i] = num / denom
    end
    return Y
end
