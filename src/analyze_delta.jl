using Statistics
using Distributions
using ProgressMeter
    conf_flag = _check_conf_flag(num_resamples, conf_level)

#=
References
----------
    [1] Borgonovo, E. (2007). "A new uncertainty importance measure."
       Reliability Engineering & System Safety, 92(6):771-784,
       doi:10.1016/j.ress.2006.04.015.

    [2] Plischke, E., E. Borgonovo, and C. L. Smith (2013). "Global
       sensitivity measures from given data." European Journal of
       Operational Research, 226(3):536-550, doi:10.1016/j.ejor.2012.11.047.
=#

"""
    analyze(data::DeltaData, model_input::AbstractArray{<:Number, S1}, model_output::AbstractArray{<:Number, S2}; num_resamples::Int = 10_000, conf_level::Number = 0.95) where S1 where S2 

Performs a Delta Moment-Independent Analysis on the `model_output` produced with 
the problem defined by the information in `data` and `model_input` and returns
a dictionary of results with the sensitivity indices and respective confidence 
intervals for each of the parameters defined using the `num_resamples` and `conf_level`
keyword args. 
"""
function analyze(data::DeltaData, model_input::AbstractArray{<:Number, S1}, model_output::AbstractArray{<:Number, S2}; num_resamples::Int = 10_000, conf_level::Number = 0.95) where S1 where S2 
    
    # this method requires a confidence interval and num_resamples, so we do not
    # currently allow these to be set to Nothing as we did with the Sobol method ...
    # if a user wants to limit run time for intial tries they can set a very low
    # number of resamples
    # conf_flag = _check_conf_flag(num_resamples, conf_level)

    # define constants
    calc_second_order = data.calc_second_order 
    D = length(data.params) # number of uncertain parameters in problem
    N = data.N # number of samples

    # values for CI calculations
    if conf_flag
        r = rand(1:N, N, num_resamples)
        Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2) # calculate z* for CI
    end

    M = min(ceil(N ^ (2 / (7 + tanh((1500 - N) / 500)))), 48)
    m = LinRange(0, N, M + 1)
    model_output_grid = LinRange(min(model_output), max(model_output), 100)

    # preallocate arrays 
    delta = Array{Float64}(undef, D)
    firstorder = Array{Float64}(undef, D)
    firstorder_conf = Array{Float64}(undef, D)
    delta_conf = Array{Float64}(undef, D)

    for i in 1:D
        delta[i], delta_conf = bias_reduced_delta(model_output, model_output_grid, model_input[:,i], m, num_resamples, conf_level)
        firstorder[i] = sobol_first(model_output, model_input[:,i], m)
        firstorder_conf[i] = sobol_first_conf(model_output, model_input[:,i], m, num_resamples, conf_level)
    end

    results = Dict(
        :delta              => delta,
        :delta_conf         => delta_conf,
        :firstorder         => firstorder,
        :firstorder_conf    => firstorder_conf
    )
    return results
end

"""
    calc_delta(model_output::AbstractArray{<:Number, S2}, model_output_grid::LinRange, 
        model_input::AbstractArray{<:Number, S1}, m::LinRange) where S1 where S2

Plischke et al. 2013 estimator (eqn 26) for d_hat
"""
function calc_delta(model_output::AbstractArray{<:Number, S2}, model_output_grid::LinRange, model_input::AbstractArray{<:Number, S1}, m::LinRange) where S1 where S2
    # TODO
end

# def calc_delta(Y, Ygrid, X, m):
#     N = len(Y)
#     fy = gaussian_kde(Y, bw_method='silverman')(Ygrid)
#     xr = rankdata(X, method='ordinal')

#     d_hat = 0
#     for j in range(len(m) - 1):
#         ix = np.where((xr > m[j]) & (xr <= m[j + 1]))[0]
#         nm = len(ix)
#         fyc = gaussian_kde(Y[ix], bw_method='silverman')(Ygrid)
#         d_hat += (nm / (2 * N)) * np.trapz(np.abs(fy - fyc), Ygrid)

#     return d_hat

"""
    bias_reduced_delta(model_output::AbstractArray{<:Number, S2}, model_output_grid::LinRange, 
        model_input::AbstractArray{<:Number, S1}, m::LinRange, num_resamples::Int, 
        conf_level::Number) where S1 where S2

Plischke et al. 2013 bias reduction technique (eqn 30)
"""
function bias_reduced_delta(model_output::AbstractArray{<:Number, S2}, model_output_grid::LinRange, model_input::AbstractArray{<:Number, S1}, m::LinRange, num_resamples::Int, conf_level::Number) where S1 where S2
    # TODO
end

# def bias_reduced_delta(Y, Ygrid, X, m, num_resamples, conf_level):
#     d = np.zeros(num_resamples)
#     d_hat = calc_delta(Y, Ygrid, X, m)

#     for i in range(num_resamples):
#         r = np.random.randint(len(Y), size=len(Y))
#         d[i] = calc_delta(Y[r], Ygrid, X[r], m)

#     d = 2 * d_hat - d
#     return (d.mean(), norm.ppf(0.5 + conf_level / 2) * d.std(ddof=1))

"""
    sobol_first(model_output::AbstractArray{<:Number, S2}, model_input::AbstractArray{<:Number, S1}, 
        m::LinRange) where S1 where S2

Definition ...
"""
function sobol_first(model_output::AbstractArray{<:Number, S2}, model_input::AbstractArray{<:Number, S1}, m::LinRange) where S1 where S2
    # TODO
end

# def sobol_first(Y, X, m):
#     xr = rankdata(X, method='ordinal')
#     Vi = 0
#     N = len(Y)
#     for j in range(len(m) - 1):
#         ix = np.where((xr > m[j]) & (xr <= m[j + 1]))[0]
#         nm = len(ix)
#         Vi += (nm / N) * (Y[ix].mean() - Y.mean()) ** 2
#     return Vi / np.var(Y)

"""
    sobol_first_conf(model_output::AbstractArray{<:Number, S2}, model_input::AbstractArray{<:Number, S1}, 
        m::LinRange, num_resamples::Int, conf_level::Number) where S1 where S2

Definition ...
"""

function sobol_first_conf(model_output::AbstractArray{<:Number, S2}, model_input::AbstractArray{<:Number, S1}, m::LinRange, num_resamples::Int, conf_level::Number) where S1 where S2
    # TODO
end

# def sobol_first_conf(Y, X, m, num_resamples, conf_level):
#     s = np.zeros(num_resamples)

#     for i in range(num_resamples):
#         r = np.random.randint(len(Y), size=len(Y))
#         s[i] = sobol_first(Y[r], X[r], m)

#     return norm.ppf(0.5 + conf_level / 2) * s.std(ddof=1)