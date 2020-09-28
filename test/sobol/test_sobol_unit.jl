using Test
using Distributions
using DataStructures

import GlobalSensitivityAnalysis: scale_samples!, ishigami

##
## 1. utils
##

# SobolData
N = 100
calc_second_order = false
parameters = OrderedDict(
    :param1 => Normal(0,1),
    :param2 => LogNormal(5, 20),
    :param3 => TriangularDist(0, 4, 1)
)

data1 = SobolData(params = parameters, calc_second_order = calc_second_order, N = N)
@test data1.params[:param1] == parameters[:param1]
@test data1.params[:param2] == parameters[:param2]
@test data1.params[:param2] == parameters[:param2]
@test data1.params[:param3] == parameters[:param3]
@test data1.N == N
@test data1.calc_second_order == calc_second_order

data2 = SobolData()
@test data2.params == nothing
@test data2.N == 1000
@test data2.calc_second_order == true

data3 = SobolData(params = parameters)
data4 = SobolData(calc_second_order = false)
data5 = SobolData(N = 100)

##
## 2. Sample Sobol
##

samples = sample(data1)
@test size(samples, 2) == length(data1.params)
if data1.calc_second_order
    @test size(samples,1) == data1.N * (2 * length(data1.params) + 2)
else
    @test size(samples,1) == data1.N * (length(data1.params) + 2)
end

samples3 = sample(data3)
@test size(samples3, 2) == length(data3.params)
if data3.calc_second_order
    @test size(samples3,1) == data3.N * (2 * length(data3.params) + 2)
else
    @test size(samples3,1) == data3.N * (length(data3.params) + 2)
end

@test_throws ErrorException sample(data2) # params are nothing
@test_throws ErrorException sample(data4) # params are nothing
@test_throws ErrorException sample(data5) # params are nothing

##
## 3a. Analyze Sobol
##

Y1 = ishigami(samples)
results = analyze(data1, Y1)
for Si in results[:firstorder]
    @test Si <= 1
end
for CI in results[:firstorder_conf]
    @test CI > 0 
end
for St in results[:totalorder]
    @test St <= 1
end
for CI in results[:totalorder_conf]
    @test CI > 0 
end
@test sum(results[:totalorder]) > sum(results[:firstorder])

Y3 = ishigami(samples3)
results = analyze(data3, Y3)
for Si in results[:firstorder]
    @test Si <= 1
end
for CI in results[:firstorder_conf]
    @test ismissing(CI) || CI > 0
end
for S2i in results[:secondorder]
    @test ismissing(S2i) || S2i <= 1
end
for CI in results[:secondorder_conf]
    @test ismissing(CI) || CI > 0
end
for St in results[:totalorder]
    @test St <= 1
end
for CI in results[:totalorder_conf]
    @test ismissing(CI) || CI > 0
end
@test sum(results[:totalorder]) > sum(results[:firstorder])

##
## 3b. Analyze Sobol Optional Keyword Args
##

data = SobolData(
    params = OrderedDict(:x1 => Normal(1, 0.2),
        :x2 => Uniform(0.75, 1.25),
        :x3 => LogNormal(0, 0.5)),
    N = 1000
)
samples = sample(data)
Y = ishigami(samples)
results = analyze(data, Y)

results = analyze(data, samples, Y; progress_meter = false) # no progress bar should show

@test length(analyze(data, samples, Y; N_override = 10)) == 4 
results_override = analyze(data, samples, Y, N_override = data.N)
results_original = analyze(data, samples, Y)
@test results_override[:firstorder] == results_original[:firstorder]
@test results_override[:totalorder] == results_original[:totalorder] 
@test_throws ErrorException analyze(data1, samples, Y1; N_override = data.N + 1) # N_override > N
