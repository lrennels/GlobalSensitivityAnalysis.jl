using Test
using Distributions
using DataStructures
using GlobalSensitivityAnalysis

import GlobalSensitivityAnalysis: ishigami, sample

##
## 1. Utils
##

# DeltaData
N = 100
parameters = OrderedDict(
    :param1 => Normal(0,1),
    :param2 => LogNormal(5, 20),
    :param3 => TriangularDist(0, 4, 1)
)

data1 = DeltaData(params = parameters, N = N)
@test data1.params[:param1] == parameters[:param1]
@test data1.params[:param2] == parameters[:param2]
@test data1.params[:param2] == parameters[:param2]
@test data1.params[:param3] == parameters[:param3]
@test data1.N == N

data2 = DeltaData()
@test data2.params == nothing
@test data2.N == 1000

data3 = DeltaData(params = parameters)
data4 = DeltaData(N = 100)

##
## 2. Sample Delta
##

samples = sample(data1)
@test size(samples, 2) == length(data1.params)
@test size(samples, 1) == data1.N

samples3 = sample(data3)
@test size(samples3, 2) == length(data3.params)
@test size(samples3, 1) == data3.N 

@test_throws ErrorException sample(data2) # params are nothing
@test_throws ErrorException sample(data4) # params are nothing

# ##
# ## 3a. Analyze Delta - Indices (TODO is there anything else we should test indices for?)
# ##

Y1 = ishigami(samples)
results = analyze(data1, samples, Y1)
for Si in results[:firstorder]
    @test Si <= 1
end
for Si in results[:delta]
    @test Si <= 1
end
for CI in results[:firstorder_conf]
    @test CI > 0 
end
for CI in results[:delta_conf]
    @test CI > 0 
end

Y3 = ishigami(samples3)
results = analyze(data3, samples3, Y3)
for Si in results[:firstorder]
    @test Si <= 1
end
for Si in results[:delta]
    @test Si <= 1
end
for CI in results[:firstorder_conf]
    @test CI > 0 
end
for CI in results[:delta_conf]
    @test CI > 0 
end

# ##
# ## 3b. Analyze Delta - Optional Keyword Args
# ##

data = DeltaData(
    params = OrderedDict(:x1 => Normal(1, 0.2),
        :x2 => Uniform(0.75, 1.25),
        :x3 => LogNormal(0, 0.5)),
    N = 1000
)
samples = sample(data)
Y = ishigami(samples)
results = analyze(data, samples, Y)

results = analyze(data, samples, Y; progress_meter = false) # no progress bar should show
results = analyze(data, samples, Y; num_resamples = 100, conf_level = 0.90, progress_meter = false) # keywords

@test length(analyze(data, samples, Y; N_override = 100)) == 5 
results_override = analyze(data, samples, Y; N_override = data.N)
results_original = analyze(data, samples, Y)
@test results_override[:firstorder] == results_original[:firstorder]
@test_throws ErrorException analyze(data, samples, Y; N_override = data.N + 1) # N_override > N
