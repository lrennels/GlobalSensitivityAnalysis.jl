using Distributions
using Random

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
    N = data.N # number of samples
    d = 1.0 / N

    # preallocate
    result = zeros(N, D)
    temp = zeros(N)

    # generate samples
    for i = 1:D

        for j = 1:N
            
            # NOTE: to perfeclty copy SALib using j - 1, probably not necessary 
            low = (j-1) * d
            high = (j-1) + 1
            uniform_distrib = Uniform(low, high)

            temp[j] = rand(uniform_distrib, 1)[1]
            shuffle!(temp)
        end

        for j = 1:N
            result[j, i] = temp[j]
        end

        scale_samples!(result, [values(data.params)...]) # scale
        return result

    end

end
