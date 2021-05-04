using Distributions
using DataStructures
using DataFrames
using Test

import GlobalSensitivityAnalysis: ishigami, split_output, first_order, total_order, sample
#
# STEP 1. try running confidence intervals a few times and looking at the means etc. 
# and then compare that to doing the same for Python SALib 
#

#
# STEP 1A. Uniform Case 
#

# define problem and some constants
N = 1_000
num_resamples = 10_000
conf_level = 0.95
calc_second_order = false
iter = 10

data = SobolData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = N,
    calc_second_order = calc_second_order
)

D = length(data.params)

# run sampling since this is deterministic
samples = sample(data) # generate samples using Sobol sequence
Y = ishigami(samples) # run model (example)

# iteratively perform Sobol Analysis and save confidence intervals
firstorder_all = zeros(iter, D)
totalorder_all = zeros(iter, D)

for i = 1:iter
    results = analyze(data, Y; num_resamples = num_resamples) 
    firstorder_all[i,:] = results[:firstorder_conf]
    totalorder_all[i,:] = results[:totalorder_conf]
end

@test mean(firstorder_all, dims = 1) ≈ [0.06369217 0.05347996 0.05821841] atol = 1e-3
@test mean(totalorder_all, dims = 1) ≈ [0.08792231 0.04021016 0.02634918] atol = 1e-3

#
# STEP 1B. Non-Uniform Case 
#

# define problem and some constants
N = 1_000
num_resamples = 10_000
conf_level = 0.95
calc_second_order = false
iter = 10

data = SobolData(
    params = OrderedDict(:x1 => Normal(1, 0.2),
        :x2 => Uniform(0.75, 1.25),
        :x3 => LogNormal(0, 0.5)),
    N = N,
    calc_second_order = calc_second_order
)

D = length(data.params)

# run sampling since this is deterministic
samples = sample(data) # generate samples using Sobol sequence
Y = ishigami(samples) # run model (example)

# iteratively perform Sobol Analysis and save confidence intervals
firstorder_all = zeros(iter, D)
totalorder_all = zeros(iter, D)

for i = 1:iter
    results = analyze(data, Y; num_resamples = num_resamples) 
    firstorder_all[i,:] = results[:firstorder_conf]
    totalorder_all[i,:] = results[:totalorder_conf]
end

@test mean(firstorder_all, dims = 1) ≈ [0.00638154 0.0948373  0.67803243] atol = 1e-1
@test mean(totalorder_all, dims = 1) ≈ [0.01245481 0.09457124 0.10595523] atol = 1e-1

#
# STEP 2. Deterministic example
#

# define problem and constants
N = 5
conf_level = 0.95
calc_second_order = false

data = SobolData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = N,
    calc_second_order = calc_second_order
)
D = length(data.params)

# sampling
samples = sample(data)
Y = ishigami(samples)

# analysis
A, B, AB, BA = split_output(Matrix(Y), N, D, calc_second_order)

# confidence intervals  - constants
r = [3 2 3;
    4 1 1;
    3 2 0;
    1 2 2;
    1 1 1]
r = r .+ 1
Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2) # calculate z* for CI

# confidence intervals - calculate
firstorder_conf = Array{Float64}(undef, D)
totalorder_conf = Array{Float64}(undef, D)

for i in 1:D
    firstorder_conf[i] = Z * std(first_order(A[r], AB[r, i], B[r]))
    totalorder_conf[i] = Z * std(total_order(A[r], AB[r, i], B[r]))
end

@test firstorder_conf ≈ [1.57157685, 0.61660013, 0.81296287] atol = 1e-8
@test totalorder_conf ≈ [2.26661086, 0.36204958, 0.46523933] atol = 1e-8
