using Distributions
using Random

# Adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
# library for sensitivity analysis. Journal of Open Source Software, 2(9)

"""
    sample(data::DeltaData)

Generate a matrix containing the model inputs for Delta-method sensitivity analysis 
with the information in the `data`. The resulting matrix has `N` rows and `D` 
columns, where `D` is the number of parameters and `N` is the number of samples.
"""
function sample(data::DeltaData)

    # check for parameters
    if data.params === nothing
        error("Cannot sample to generate model inputs because there are no parameters")
    end

    # constants     
    D = length(data.params) # number of uncertain parameters in problem
    d = 1.0 / data.N

    # preallocate
    result = zeros(data.N, D)
    temp = zeros(data.N)

    # generate samples
    for i = 1:D
        for j = 1:data.N
            
            low = (j - 1) * d
            high = j * d
            uniform_distrib = Uniform(low, high)

            temp[j] = rand(uniform_distrib)
        end
        shuffle!(temp)
        result[:, i] = temp
    end

    scale_samples!(result, [values(data.params)...]) # scale
    return result

end
