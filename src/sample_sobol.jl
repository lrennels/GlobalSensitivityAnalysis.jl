using Distributions

include("utils.jl")
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
=#

"""
    sample(data::SobolPayload)

Generate a matrix containing the model inputs for Sobol sensitivity analysis with 
the information in the `params` payload. In this function we apply Saltelli's 
extension of the Sobol  sequence. Saltelli's scheme extends the Sobol sequence in 
a way to reduce the error rates in the resulting sensitivity index calculations. 
The resulting matrix has `N` * (`D` + 2) rows, where `D` is the number of parameters 
and `N` is the number of samples.
"""
# TODO - include second order effects
# TODO - include groups
function sample(data::SobolPayload)

    # set number of values to skip from the initial sequence 
    numskip = 1000

    # constants 
    D = length(data.params) # number of uncertain parameters in problem
    N = data.N # number of samples

    base_seq = sobol_sequence(N + numskip, 2 * D)
    base_seq = scale_sobol_seq(base_seq, [values(data.params)...]) #scale

    index = 1

    # create the Saltelli sequence
    saltelli_seq = Array{Float64}(undef, N * (D + 2), D)

    # The Saltelli sequence is made up of N blocks of (D + 2) rows, where each block
    # j contains a first row from A_j, a last row B_j and D middle rows that form 
    # AB_j. [Saltelli et al., 2010 see Radial Sampling]

    for i in (numskip + 1): (N + numskip)

        # copy matrix "A" (first row of each block)
        for j in 1:D
            saltelli_seq[index, j] = base_seq[i, j]
        end
        index += 1

        # for each parameter, place elements of "B" into "A" and insert those D rows
        # of "AB" (middle rows of each block)
        for k in 1:D
            for j in 1:D
                if j == k
                    saltelli_seq[index, j] = base_seq[i, j + D] # from B
                else
                    saltelli_seq[index, j] = base_seq[i, j] # from A
                end
            end
            index += 1
        end

        # copy matrix "B" (last row of each block)
        for j in 1:D
            saltelli_seq[index, j] = base_seq[i, j + D]
        end
        index += 1
    end

    return saltelli_seq
end
