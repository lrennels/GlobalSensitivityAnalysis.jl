using Distributions
using DataStructures
using DataFrames
using Test
using KernelDensity
using NumericalIntegration
using Statistics

import StatsBase: ordinalrank
import GlobalSensitivityAnalysis: calc_delta, sobol_first

#
# STEP 1. Deterministic example
#

# define problem and constants
N = 5
num_resamples = 3
conf_level = 0.95

data = DeltaData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = N
)

# sampling and eval
model_input =   [1.01124        -2.32431    1.52026
            -3.10386        -0.222252   1.98626
            2.52807         2.21401     -0.920672
            -0.000432765    0.851149    -0.559435
            -1.87759        -1.59184    -2.88272]

model_output = ishigami(convert(Matrix, model_input))

# some setup
if size(model_output, 2) != 1
    error("Model output for analyzing DeltaData has more than one col, not handled yet.")
else
    model_output = vec(model_output)
end

D = length(data.params) # number of uncertain parameters in problem
M = Int(min(ceil(N ^ (2 / (7 + tanh((1500 - N) / 500)))), 48))
m = collect(LinRange(0, N, M + 1))
model_output_grid = collect(LinRange(minimum(model_output), maximum(model_output), 100))

firstorder_conf = Array{Float64}(undef, D)
delta_conf = Array{Float64}(undef, D)

# calculate confidence intervals

for i in 1:D

    #
    # delta
    #
    
    d = zeros(num_resamples)
    d_hat = calc_delta(model_output, model_output_grid, model_input[:,i], m)
    for j = 1:num_resamples

        # r = np.random.randint(len(Y), size=len(Y))
        if j == 1
            r = [1,2,3,4,5]
        elseif j == 2
            r = [3,4,5,1,1]
        else
            r = [4,1,2,5,5]
        end

        d[j] = calc_delta(model_output[r], model_output_grid, model_input[r,i], m)
    end

    d = 2 * d_hat .- d
    Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2)
    delta_conf[i] = Z * std(d)

    #
    # first order
    #

    s = zeros(num_resamples)
    Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2)

    for j in 1:num_resamples
        # r = np.random.randint(len(Y), size=len(Y))
        if j == 1
            r = [1,2,3,4,5]
        elseif j == 2
            r = [3,4,5,1,1]
        else
            r = [4,1,2,5,5]
        end
        s[j] = sobol_first(model_output[r], model_input[r,i], m)
    end

    firstorder_conf[i] = Z * std(s)

end
