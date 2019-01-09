using Distributions

include("sobol_sequence.jl")

#=
Code adapted from: Herman, J. and Usher, W. (2017) SALib: An open-source Python 
library for sensitivity analysis. Journal of Open Source Software, 2(9)

References
----------
    [1] Sobol, I. M. (2001).  "Global sensitivity indices for nonlinear
        mathematical models and their Monte Carlo estimates."  Mathematics
        and Computers in Simulation, 55(1-3):271-280,
        doi:10.1016/S0378-4754(00)00270-6.
    [2] Saltelli, A. (2002).  "Making best use of model evaluations to
        compute sensitivity indices."  Computer Physics Communications,
        145(2):280-297, doi:10.1016/S0010-4655(02)00280-1.
    [3] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
        S. Tarantola (2010).  "Variance based sensitivity analysis of model
        output.  Design and estimator for the total sensitivity index."
        Computer Physics Communications, 181(2):259-270,
        doi:10.1016/j.cpc.2009.09.018.
    [4] Saltelli, Andrea, et al. Global sensitivity analysis: the primer. 
        John Wiley & Sons, 2008.
=#

"""
    SobolParams

Implement a `SobolParams` type containing the `names` of the uncertain parameters
of the problem, and the distributions `dists`, corresponding to each parameter.
"""
mutable struct SobolParams
    names::Array{String, 1}
    dists::Array{Distribution, 1}
end

"""
    scale_sobol_seq(sequence::Array, dists::Array{Distribution, 1})

Rescale a Sobol `sequence` of parameters from the 0-to-1 range to the defined bounds
of their distributions `dists`.  This function requires the parameters to have
uniform distributions as of now.
"""
function scale_sobol_seq(sequence::Array, dists::Array{Distribution, 1})
    for dist in dists
        if !(typeof(dist)  <: Uniform)
            error("scale_params can only scale parameters with Uniform distributions")
        end
    end
    sequence = sequence .* (repeat(scale.(dists), 2))' .+ (repeat(minimum.(dists), 2))'

    return sequence
end

"""
    saltelli_sample(params::SobolParams, N::Int, method="quasi-random")

Generate a matrix containing the model inputs for Sobol sensitivity analysis with `N` 
samples and uncertain parameters described by `params`. The method used to
generate the sequence is defined by `seq_method` and can be (1) `quasi-random`, 
currently the Sobol sequence, or (2) `random`. Regardless of the method used, we
then apply Saltelli's extension of the Sobol sequence. Saltelli's scheme extends 
the Sobol sequence in a way to reduce the error rates in the resulting sensitivity 
index calculations, which is irrelevant if the `.  

The resulting matrix has `N` * (`D` + 2) rows, where `D` is the number of parameters.  
"""
# TODO - include second order effects
# TODO - include groups
function saltelli_sample(params::SobolParams, N::Int, method="quasi-random")

    # set number of values to skip from the initial sequence (TODO - why do we need 
    # to skip values?  used so can validate against SALib but unclear why this is required)
    numskip = 1000

    # number of uncertain parameters in problem
    D = length(params.names)

    if method == "quasi-random"
        base_seq = sobol_sequence(N + numskip, 2 * D)
        base_seq = base_seq[numskip + 1:end, :]

        # scale 
        base_seq = scale_sobol_seq(base_seq, params.dists)
    else
        # TODO: how to we want to plug random in here?  Mimi hook.
    end

    # create the Saltelli sequence
    saltelli_seq = zeros(N * (D + 2), D)

    # copy matrix "A" into the beginnning of the saltelli sequence [4]
    saltelli_seq[1:N, :] = base_seq[1:N, 1:D]

    # for each parameter, place elements of "B" into "A" and insert into Saltelli
    # sequence to create an "AB" for each parameter [4]
    for i in 1:D
        row = i * N + 1
        saltelli_seq[row:(row + N - 1), :] = base_seq[:, 1:D] # AB = A
        saltelli_seq[row:(row + N - 1), i] = base_seq[:, D + i] # insert slice of B
    end

    # copy matrix "B" into the end of the saltelli sequence [4]
    saltelli_seq[end - N + 1:end, :] = base_seq[:, D+1:end]

    # scale the distributions
    return saltelli_seq
end