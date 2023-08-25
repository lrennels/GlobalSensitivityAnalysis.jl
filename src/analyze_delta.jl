import StatsBase: ordinalrank

# Adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
# library for sensitivity analysis. Journal of Open Source Software, 2(9)

# References
# ----------
#     [1] Borgonovo, E. (2007). "A new uncertainty importance measure."
#        Reliability Engineering & System Safety, 92(6):771-784,
#        doi:10.1016/j.ress.2006.04.015.
#     [2] Plischke, E., E. Borgonovo, and C. L. Smith (2013). "Global
#        sensitivity measures from given data." European Journal of
#        Operational Research, 226(3):536-550, doi:10.1016/j.ejor.2012.11.047.

"""
    analyze(data::DeltaData, 
            model_input::AbstractArray{<:Number, S1}, 
            model_output::AbstractArray{<:Number, S2}; 
            num_resamples::Int = 1_000, 
            conf_level::Number = 0.95, 
            progress_meter::Bool = true, 
            N_override::Union{Nothing, Integer}=nothing) where S1 where S2 

Performs a Delta Moment-Independent Analysis on the `model_output` produced with 
the problem defined by the information in `data` and `model_input` and returns
a dictionary of results with the sensitivity indices and respective confidence 
intervals for each of the parameters defined using the `num_resamples` and `conf_level`
keyword args.  The `progress_meter` keyword argument indicates whether a progress meter 
will be displayed and defaults to true. The `N_override` keyword argument allows users 
to override the `N` used in a specific `analyze` call to analyze just a subset 
(useful for convergence graphs).
"""
function analyze(data::DeltaData,
                model_input::AbstractArray{<:Number, S1}, 
                model_output::AbstractArray{<:Number, S2}; 
                num_resamples::Int = 1_000, 
                conf_level::Number = 0.95, 
                progress_meter::Bool = true, 
                N_override::Union{Nothing, Integer}=nothing) where S1 where S2 
    
    # This method requires a confidence interval and num_resamples, so we do not
    # currently allow these to be set to Nothing as we did with the Sobol method.
    # If a user wants to limit run time for intial tries they can set a very low
    # number of resamples

    @warn "The GlobalSensitivityAnalysis delta method is a beta-release undergoing further testing, please use the method with this knowledge and notify the developers with any issues or questions."

    if size(model_output, 2) != 1
        error("Model output for analyzing DeltaData has more than one column, this is not yet supported.")
    else
        model_output = vec(model_output) # convert to a Vector
    end

    # define constants
    D = length(data.params) # number of uncertain parameters in problem

    # deal with overriding N
    if isnothing(N_override)
        N = data.N # number of samples
    else
        N_override > data.N ? error("N_override ($N_override) cannot be greater than original N used in sampling ($(data.N))") : nothing 
        N = N_override # number of samples

        # reduce the output to just what should be considered for this N
        model_output = model_output[1:N]
        model_input = model_input[1:N,:]
    end

    exp =  (2 / (7 + tanh((1500 - N) / 500)))
    M = Int(min(ceil(N ^ exp), 48))
    m = collect(LinRange(0, N, M + 1))
    model_output_grid = collect(LinRange(minimum(model_output), maximum(model_output), 100)) # note that SALIb defaults the YGridlength to 100, while GlobalSensitivity.jl defaults to 2048

    # preallocate arrays 
    delta_biased = Array{Float64}(undef, D)
    delta = Array{Float64}(undef, D)
    delta_conf = Array{Float64}(undef, D)
    firstorder = Array{Float64}(undef, D)
    firstorder_conf = Array{Float64}(undef, D)

    # set up progress meter
    counter = 0
    progress_meter ? p = Progress(D, counter, "Calculating indices for $D parameters ...") : nothing

    for i in 1:D
        # increment progress meter
        counter += 1
        progress_meter ? ProgressMeter.update!(p, counter) : nothing  
        delta_biased[i], delta[i], delta_conf[i] = bias_reduced_delta(model_output, model_output_grid, model_input[:,i], m, num_resamples, conf_level)
        firstorder[i] = sobol_first(model_output, model_input[:,i], m)
        firstorder_conf[i] = sobol_first_conf(model_output, model_input[:,i], m, num_resamples, conf_level)
    end

    results = Dict(
        :delta_biased       => delta_biased,
        :delta              => delta, 
        :delta_conf         => delta_conf,
        :firstorder         => firstorder,
        :firstorder_conf    => firstorder_conf
    )
    return results
end

"""
    calc_delta(model_output::AbstractArray{<:Number, S1}, model_output_grid::AbstractArray{<:Number, 1}, 
                model_input::AbstractArray{<:Number, S2}, m::AbstractArray{<:Number, 1}) where S1 where S2

Plischke et al. 2013 estimator (eqn 26) for d_hat
"""
function calc_delta(model_output::AbstractArray{<:Number, S1}, model_output_grid::AbstractArray{<:Number, 1}, 
                    model_input::AbstractArray{<:Number, S2}, m::AbstractArray{<:Number, 1}) where S1 where S2
    N = length(model_output)
    k = KernelDensity.kde(model_output) # defaults are kernel = normal and bandwidth = Silverman which match SALib specification
    fy = Distributions.pdf(k, model_output_grid) # eq 23.1
    model_input_ranks = ordinalrank(model_input)

    d_hat = 0
    for j = 1:length(m) - 1
        ix = findall((model_input_ranks .> m[j]) .& (model_input_ranks .<= m[j + 1]))
        nm = length(ix)

        # get the separation between the total y pdf and the condition y pdf induced by this class

        # Estimated conditional distribution of y (using normal kernel)
        k = KernelDensity.kde(model_output[ix]) # defaults are kernel = normal and bandwidth = Silverman which match SALib
        fyc = Distributions.pdf(k, model_output_grid) # eq 23.2
        pdf_diff = abs.(fy - fyc) # eq 24

        # Use trapezoidal rule to estimate the difference between the curves.
        class_separation = Trapz.trapz(model_output_grid, pdf_diff) # eq 25

        # Increment Estimator
        d_hat += (nm / (2 * N)) * class_separation # eq 26

    end
    return d_hat
    
end

"""
    bias_reduced_delta(model_output::AbstractArray{<:Number, S1}, model_output_grid::AbstractArray{<:Number, 1}, 
                        model_input::AbstractArray{<:Number, S2}, m::AbstractArray{<:Number, 1}, num_resamples::Int, 
                        conf_level::Number) where S1 where S2

Plischke et al. 2013 bias reduction technique (eqn 30)
"""
function bias_reduced_delta(model_output::AbstractArray{<:Number, S1}, model_output_grid::AbstractArray{<:Number, 1}, 
                            model_input::AbstractArray{<:Number, S2}, m::AbstractArray{<:Number, 1}, num_resamples::Int,
                            conf_level::Number) where S1 where S2

    d = zeros(num_resamples)
    d_hat = calc_delta(model_output, model_output_grid, model_input, m)

    for i = 1:num_resamples
        r = rand(1:length(model_output), length(model_output))
        d[i] = calc_delta(model_output[r], model_output_grid, model_input[r], m)
    end

    d = 2 * d_hat .- d
    Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2)

    return (d_hat, mean(d), Z * std(d))
end

"""
    sobol_first(model_output::AbstractArray{<:Number, S1}, model_input::AbstractArray{<:Number, S2}, 
                m::AbstractArray{<:Number, 1}) where S1 where S2
"""
function sobol_first(model_output::AbstractArray{<:Number, S1}, model_input::AbstractArray{<:Number, S2}, 
                    m::AbstractArray{<:Number, 1}) where S1 where S2
    model_input_ranks = ordinalrank(model_input)
    Vi = 0
    N = length(model_output)
    for j in 1:length(m) - 1
        ix = findall((model_input_ranks .> m[j]) .& (model_input_ranks .<= m[j + 1]))
        nm = length(ix)
        Vi += (nm / N) * ((mean(model_output[ix]) - mean(model_output))^2)
    end

    return Vi / var(model_output)

end

"""
    sobol_first_conf(model_output::AbstractArray{<:Number, S2}, model_input::AbstractArray{<:Number, S1}, 
        m::AbstractArray{<:Number, 1}, num_resamples::Int, conf_level::Number) where S1 where S2
"""
function sobol_first_conf(model_output::AbstractArray{<:Number, S2}, model_input::AbstractArray{<:Number, S1},
                        m::AbstractArray{<:Number, 1}, num_resamples::Int, conf_level::Number) where S1 where S2
                        
    s = zeros(num_resamples)
    Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2)

    for i in 1:num_resamples
        r = rand(1:length(model_output), length(model_output))
        s[i] = sobol_first(model_output[r], model_input[r], m)
    end

    return Z * std(s)
end
