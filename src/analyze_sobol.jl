# Adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
# library for sensitivity analysis. Journal of Open Source Software, 2(9)

# References
#
#     [1] Sobol, I. M. (2001).  "Global sensitivity indices for nonlinear
#         mathematical models and their Monte Carlo estimates."  Mathematics
#         and Computers in Simulation, 55(1-3):271-280,
#         doi:10.1016/S0378-4754(00)00270-6.
#     [2] Saltelli, A. (2002).  "Making best use of model evaluations to
#         compute sensitivity indices."  Computer Physics Communications,
#         145(2):280-297, doi:10.1016/S0010-4655(02)00280-1.
#     [3] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
#         S. Tarantola (2010).  "Variance based sensitivity analysis of model
#         output.  Design and estimator for the total sensitivity index."
#         Computer Physics Communications, 181(2):259-270,
#         doi:10.1016/j.cpc.2009.09.018.

"""
    analyze(data::SobolData, 
            model_output::AbstractArray{<:Number, S}; 
            num_resamples::Union{Nothing, Int} = 1_000, 
            conf_level::Union{Nothing, Number} = 0.95, 
            progress_meter::Bool = true,
            N_override::Union{Nothing, Integer}=nothing) where S

Performs a Sobol Analysis on the `model_output` produced with the problem 
defined by the information in `data` and returns the a dictionary of results
with the sensitivity indices and respective confidence intervals for each of the
parameters defined using the `num_resamples` and `conf_level` keyword args. If these
are nothing than no confidence intervals will be calculated. The `progress_meter`
keyword argument indicates whether a progress meter will be displayed and defaults
to true. The `N_override` keyword argument allows users to override the `N` used in
a specific `analyze` call to analyze just a subset (useful for convergence graphs).
"""
function analyze(data::SobolData, 
                model_output::AbstractArray{<:Number, S}; 
                num_resamples::Union{Nothing, Int} = 1_000, 
                conf_level::Union{Nothing, Number} = 0.95, 
                progress_meter::Bool = true, 
                N_override::Union{Nothing, Integer}=nothing) where S

    # handle confidence interval flag
    if isnothing(conf_level)
        conf_flag = false
    elseif !isnothing(num_resamples)
        conf_flag = true
    else
        error("A non-nothing confidence level ($conf_level) requires a specified number of resamples.")
    end

    # define constants
    D = length(data.params) # number of uncertain parameters in problem

    # deal with overriding N
    if N_override === nothing
        N = data.N # number of samples
    else
        N_override > data.N ? error("N_override ($N_override) cannot be greater than original N used in sampling ($(data.N))") : nothing 
        N = N_override # number of samples

        # reduce the output to just what should be considered for this N
        if data.calc_second_order
            lastrow = N * ((2 * D) + 2)
        else
            lastrow = N * (D + 2)
        end
        model_output = model_output[1:lastrow]
    end

    # values for CI calculations
    if conf_flag
        r = rand(1:N, N, num_resamples)
        Z = quantile(Normal(0.0, 1.0),1 - (1 - conf_level)/2) # calculate z* for CI
    end
    
    # normalize model output
    model_output = (model_output .- mean(model_output)) ./ std(model_output)

    # separate the model_output into results from matrices "A". "B" and "AB" 
    A, B, AB, BA = split_output(model_output, N, D, data.calc_second_order)

    # preallocate arrays for indices
    firstorder = Array{Float64}(undef, D)
    totalorder = Array{Float64}(undef, D)
    data.calc_second_order ? secondorder =  fill!(Array{Union{Float64, Missing}}(undef, D, D), missing) : nothing

    # preallocate arrays for confidence intervals
    if conf_flag
        firstorder_conf = Array{Float64}(undef, D)
        totalorder_conf = Array{Float64}(undef, D)
        data.calc_second_order ? secondorder_conf =  fill!(Array{Union{Float64, Missing}}(undef, D, D), missing) : nothing
    end

    # set up progress meter
    counter = 0
    progress_meter ? p = Progress(D, counter, "Calculating indices for $D parameters ...") : nothing

    for i in 1:D

        # increment progress meter
        counter += 1
        progress_meter ? ProgressMeter.update!(p, counter) : nothing  
        
        # first order and total order indices
        firstorder[i] = first_order(A, AB[:, i], B)[1] # array to scalar with [1]
        totalorder[i] = total_order(A, AB[:, i], B)[1] # array to scalar with [1]

        # first order and total order indice confidence intervals
        conf_flag ? firstorder_conf[i] = Z * std(first_order(A[r], AB[r, i], B[r])) : nothing
        conf_flag ? totalorder_conf[i] = Z * std(total_order(A[r], AB[r, i], B[r])) : nothing

        # second order indices
        if data.calc_second_order
            for j in (i+1):D
                secondorder[i, j] = second_order(A, AB[:, i], AB[:, j], BA[:, i], B)[1] # array to scalar with [1]
                conf_flag ? secondorder_conf[i,j] = Z * std(skipmissing(second_order(A[r], AB[r, i], AB[r, j], BA[r, i], B[r]))) : nothing
            end
        end
    end

    if data.calc_second_order
        if conf_flag
            results = Dict(
                :firstorder         => firstorder,
                :firstorder_conf    => firstorder_conf,
                :totalorder         => totalorder,
                :totalorder_conf    => totalorder_conf,
                :secondorder        => secondorder,
                :secondorder_conf   => secondorder_conf
            )
        else
            results = Dict(
                :firstorder         => firstorder,
                :totalorder         => totalorder,
                :secondorder        => secondorder,
            )
        end
    else 
        if conf_flag
            results = Dict(
                :firstorder         => firstorder,
                :firstorder_conf    => firstorder_conf,
                :totalorder         => totalorder,
                :totalorder_conf    => totalorder_conf
            )
        else
            results = Dict(
                :firstorder         => firstorder,
                :totalorder         => totalorder,
            )
        end
    end
    return results
end

"""
    first_order(A::AbstractArray{<:Number, N}, AB::AbstractArray{<:Number, N}, B::AbstractArray{<:Number, N})

Calculate the first order sensitivity indices for model outputs given model outputs
separated out into `A`, `AB`, and `A` and normalize by the variance of `[A B]`. [Saltelli et al., 
2010 Table 2 eq (b)]
"""
function first_order(A::AbstractArray{<:Number, N}, AB::AbstractArray{<:Number, N}, B::AbstractArray{<:Number, N}) where N
    return (mean(B .* (AB .- A), dims = 1) ./ var(vcat(A, B), dims = 1, corrected = false))
end

"""
    second_order(A::AbstractArray{<:Number, N}, ABi::AbstractArray{<:Number, N}, ABj::AbstractArray{<:Number, N}, BAi::AbstractArray{<:Number, N}, B::AbstractArray{<:Number, N}) where N

Calculate the second order sensitivity index between two parameters for model outputs 
given model outputs separated out into `A`, `AB`, `BA`, and `B` and normalize by 
the variance of `[A B]`. [Saltelli et al. , 2002]
"""
function second_order(A::AbstractArray{<:Number, N}, ABi::AbstractArray{<:Number, N}, ABj::AbstractArray{<:Number, N}, BAi::AbstractArray{<:Number, N}, B::AbstractArray{<:Number, N}) where N

    Vj = (mean(BAi .* ABj .- A .* B, dims = 1) ./ var(vcat(A, B), dims = 1, corrected = false))
    Si = first_order(A, ABi, B)
    Sj = first_order(A, ABj, B)

    return Vj .- Si .- Sj
end

"""
    total_order(A::AbstractArray{<:Number, N}, AB::AbstractArray{<:Number, N}, B::AbstractArray{<:Number, N})

Calculate the total order sensitivity indices for model outputs given model outputs
separated out into `A`, `AB`, and `A` and normalize by the variance of `[A B]`. [Saltelli et al., 
2010 Table 2 eq (f)].
"""
function total_order(A::AbstractArray{<:Number, N}, AB::AbstractArray{<:Number, N}, B::AbstractArray{<:Number, N}) where N
    return (0.5 * mean((A .- AB).^2, dims = 1) ./ var(vcat(A, B), dims = 1, corrected = false))
end

"""
    split_output(model_output::AbstractArray{<:Number, S}, N, D, calc_second_order)

Separate the `model_outputs` into matrices "A", "B", "AB", and "BA" for calculation 
of sensitvity indices and return those four matrices. If `calc_second_order` is 
`False`, `BA` will be `nothing`.
"""
function split_output(model_output::AbstractArray{<:Number, S}, N, D, calc_second_order::Bool) where S

    if calc_second_order
        stepsize = 2 * D + 2 
    else
        stepsize = D + 2
    end

    A = model_output[1:stepsize:end]
    B = model_output[stepsize:stepsize:end]
    
    #preallocate
    AB = Array{Float64}(undef, N, D)
    if calc_second_order
        BA = Array{Float64}(undef, N, D)
    else
        BA = nothing
    end

    for i in 1:D
        AB[:, i] = model_output[i+1:stepsize:end, :] 
        if calc_second_order
            BA[:, i] = model_output[i + D + 1:stepsize:end, :]
        end
    end

    return A, B, AB, BA
end
