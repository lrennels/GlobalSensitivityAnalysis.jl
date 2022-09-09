# file compied from https://github.com/davidanthoff/WRS.jl/tree/master/src/WRS.jl
# with permission from Professor David Anthoff, license terms as follows:

# The WRS.jl package is licensed under the MIT "Expat" License:

# Copyright (c) 2018-2019: David Anthoff and Cora Kingdon.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

# method used in Moore, F. C., Rising, J., Lollo, N., Springer, C., Vasquez, V., 
# Dolginow, A., ... & Anthoff, D. (2018). Mimi-PAGE, an open-source implementation 
# of the PAGE09 integrated assessment model. Scientific data, 5(1), 1-8.


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
