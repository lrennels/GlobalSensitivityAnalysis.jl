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
    names::AbstractArray{String, 1}
    dists::AbstractArray{Distribution, 1}
end

"""
    scale_sobol_seq(sequence::AbstractArray{<:Number, 2}, dists::AbstractArray{Distribution, 1})

Rescale a Sobol `sequence` of parameters from the 0-to-1 range to the defined bounds
of their distributions `dists`.  This function requires the parameters to have
uniform distributions as of now.
"""
function scale_sobol_seq(sequence::AbstractArray{<:Number, 2}, dists::AbstractArray{Distribution, 1})
    for dist in dists
        if !(typeof(dist)  <: Uniform)
            error("scale_params can only scale parameters with Uniform distributions")
        end
    end
    sequence = sequence .* (repeat(scale.(dists), 2))' .+ (repeat(minimum.(dists), 2))'

    return sequence
end

"""
    saltelli_sample(params::SobolParams, N::Int)

Generate a matrix containing the model inputs for Sobol sensitivity analysis with `N` 
samples and uncertain parameters described by `params`. We then apply Saltelli's 
extension of the Sobol  sequence. Saltelli's scheme extends the Sobol sequence in 
a way to reduce the error rates in the resulting sensitivity index calculations. 
The resulting matrix has `N` * (`D` + 2) rows, where `D` is the number of parameters. 
"""
# TODO - include second order effects
# TODO - include groups
function saltelli_sample(params::SobolParams, N::Int)

    # set number of values to skip from the initial sequence 
    numskip = 1000

    # number of uncertain parameters in problem
    D = length(params.names)

    base_seq = sobol_sequence(N + numskip, 2 * D)
    base_seq = scale_sobol_seq(base_seq, params.dists) #scale

    index = 1

    # create the Saltelli sequence
    saltelli_seq = zeros(N * (D + 2), D)

    for i in (numskip + 1): (N + numskip)

        # copy matrix "A" into the beginnning of the saltelli sequence [4]
        for j in 1:D
            saltelli_seq[index, j] = base_seq[i, j]
        end
        index += 1

        # for each parameter, place elements of "B" into "A" and insert into Saltelli
        # sequence to create an "AB" for each parameter [4]
        for k in 1:D
            for j in 1:D
                if j == k
                    saltelli_seq[index, j] = base_seq[i, j + D]
                else
                    saltelli_seq[index, j] = base_seq[i, j]
                end
            end
            index += 1
        end

        # copy matrix "B" into the the saltelli sequence [4]
        for j in 1:D
            saltelli_seq[index, j] = base_seq[i, j + D]
        end
        index += 1
    end

    # return
    return saltelli_seq
end