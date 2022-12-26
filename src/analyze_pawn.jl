using Logging
using Statistics, Distributions, HypothesisTests
using ProgressMeter


"""
    ks_statistic(ks)

Calculate the Kolmogorov-Smirnov test statistic.
"""
function ks_statistic(ks::HypothesisTests.KSTest)::Float64
    n = (ks.n_x * ks.n_y) / (ks.n_x + ks.n_y)

    return sqrt(n) * ks.δ
end


"""
    analyze(
        data::PAWNData,
        model_input::AbstractArray{<:Number, S1},
        model_output::AbstractArray{<:Number, S2};
        progress_meter::Bool = true
    )::Dict{Symbol, Vector{Float64}} where S1 where S2

Calculates the PAWN sensitivity index.

# Arguments
- data
- model_input
- model_output
- progress_meter

# Returns
Dict, of min, mean, median, max, std, and cv summary statistics.

# References
1. Pianosi, F., Wagener, T., 2018.
   Distribution-based sensitivity analysis from a generic input-output sample.
   Environmental Modelling & Software 108, 197-207.
   https://doi.org/10.1016/j.envsoft.2018.07.019

2. Baroni, G., Francke, T., 2020.
   GSA-cvd
   Combining variance- and distribution-based global sensitivity analysis
   https://github.com/baronig/GSA-cvd
"""
function analyze(
    data::PAWNData,
    model_input::AbstractArray{<:Number, S1},
    model_output::AbstractArray{<:Number, S2};
    progress_meter::Bool = true
)::Dict{Symbol, Vector{Float64}} where S1 where S2

    S = data.S

    N, D = size(model_input)
    step = 1 / S

    X_di = zeros(N)
    X_q = zeros(S + 1)
    pawn_t = zeros(S + 1, D)
    t_res = zeros(D, 6)

    seq = 0:step:1

    counter = 0
    progress_meter ? p = Progress(D, counter, "Calculating indices for $D parameters ...") : nothing

    # Hide warnings from HypothesisTests about ties
    with_logger(NullLogger()) do
        for d_i in 1:D
            counter += 1
            progress_meter ? ProgressMeter.update!(p, counter) : nothing  

            X_di .= model_input[:, d_i]
            X_q .= quantile(X_di, seq)
            for s in 1:S
                Y_sel = model_output[(X_di.>=X_q[s]).&(X_di.<X_q[s+1])]
                if length(Y_sel) == 0
                    pawn_t[s, d_i] = 0.0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, model_output))
            end

            p_ind = pawn_t[:, d_i]
            p_mean = mean(p_ind)
            p_sdv = std(p_ind)
            p_cv = p_sdv ./ p_mean
            t_res[d_i, :] .= (minimum(p_ind), p_mean, median(p_ind), maximum(p_ind), p_sdv, p_cv)
        end
    end

    replace!(t_res, NaN => 0.0, Inf => 0.0)

    results::Dict{Symbol, Vector{Float64}} = Dict(
        :min => t_res[:, 1],
        :mean => t_res[:, 2],
        :median => t_res[:, 3],
        :max => t_res[:, 4],
        :std => t_res[:, 5],
        :cv => t_res[:, 6],
    )

    return results
end
