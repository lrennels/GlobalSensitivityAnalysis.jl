# file compied from https://github.com/davidanthoff/WRS.jl/tree/master/src/WRS.jl

module WRS

using Distributions
using Distributed

export hd, pb2gen

# This implements (3.16) from (2017) p. 71-72
function hd(X, q=0.5, issorted=false)
    n = length(X)

    a = (n+1)*q
    b = (n+1)*(1-q)

    β = Beta(a, b)

    sorted_X = issorted ? X : sort(X)
    
    θ = 0.
    for i=1:n
        w = cdf(β, i/n) - cdf(β, (i-1)/n)
        θ += w*sorted_X[i]
    end

    return θ
end

function pb2gen(x, y; quantiles=[0.05, 0.25, 0.5, 0.75, 0.95], alpha=0.05, nboot=2000, estimator=hd)
    bootstrapped_diff_est = pmap(1:nboot) do i
        sampled_values_x = sample(x, length(x), replace=true)
        sampled_values_y = sample(y, length(y), replace=true)

        sort!(sampled_values_x)
        sort!(sampled_values_y)

        return [estimator(sampled_values_x, q, true) - estimator(sampled_values_y, q, true) for q in quantiles]
    end
    bootstrapped_diff_estimates = [[bootstrapped_diff_est[r_i][q_i] for r_i in 1:nboot] for q_i in 1:length(quantiles)]

    return map(enumerate(bootstrapped_diff_estimates)) do (q_i, bootstrapped_diff_est)
        sort!(bootstrapped_diff_est)

        low = round(Int, (alpha/2)*nboot)+1
        up = nboot-low

        q = quantiles[q_i]

        est_x = estimator(x, q)
        est_y = estimator(y, q)
        est_diff = est_x - est_y
        ci_low = bootstrapped_diff_est[low]
        ci_up = bootstrapped_diff_est[up]

        A = count(i->i<0,bootstrapped_diff_est)
        C = count(i->i==0,bootstrapped_diff_est)    
        p_hat_star = A/nboot+0.5*C/nboot
        pvalue = 2*(min(p_hat_star,1-p_hat_star))

        se = var(bootstrapped_diff_est)

        return (quantile=q, est_x=est_x, est_y=est_y, est_diff=est_diff, ci_low=ci_low, ci_up=ci_up, pvalue=pvalue, se=se, signif=pvalue < alpha)
    end
end

end # module
