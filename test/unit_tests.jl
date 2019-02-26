using Test
using Distributions

include("../src/utils.jl")
include("../src/sample_sobol.jl")
include("../src/test_functions/ishigami.jl")
include("../src/analyze_sobol.jl")

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

# scale_sobol_seq!
seq = rand(100, 8)
original_seq = copy(seq)
dists = [Normal(1, 0.2), Uniform(0.75, 1.25), LogNormal(20,4), TriangularDist(0, 4, 1)]
scale_sobol_seq!(seq, dists)
@test size(seq) == size(original_seq)
@test seq != original_seq

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

@test_throws MethodError samples2 = sample(data2) # params are nothing
@test_throws MethodError samples4 = sample(data4) # params are nothing
@test_throws MethodError samples5 = sample(data5) # params are nothing

##
## 4. Analyze Sobol
##

Y1 = ishigami(samples)
results1 = analyze(data1, Y1)
for Si in results1[:firstorder]
    @test Si <= 1
end
for St in results1[:totalorder]
    @test St <= 1
end
@test sum(results1[:totalorder]) > sum(results1[:firstorder])

Y3 = ishigami(samples3)
results3 = analyze(data3, Y3)
for Si in results3[:firstorder]
    @test Si <= 1
end
for S2i in results3[:secondorder]
    @test ismissing(S2i) || S2i <= 1
end
for St in results3[:totalorder]
    @test St <= 1
end
@test sum(results3[:totalorder]) > sum(results3[:firstorder])
