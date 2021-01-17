using VegaLite
using GlobalSensitivityAnalysis
using Distributions
using DataFrames
using CSVFiles
using DataStructures

import GlobalSensitivityAnalysis: ishigami

# Goal: Replicate Fig 3. in Plischke et al. 2013
# the SciML comparison is here in Section 4: https://github.com/danielalfonsetti/18337FinalProject/blob/master/18337final_project_daniel_alfonsetti.pdf

output_dir = "/Users/lisarennels/.julia/dev/GlobalSensitivityAnalysis/output/LHS_QuantileTesting/Plischke2013_Fig3"
mkdir(output_dir)

sample_sizes = 2 .^ [i for i in 9:14]
fig3_biased_deltas = zeros(length(sample_sizes))
fig3_deltas = zeros(length(sample_sizes))
fig3_deltas_CI = zeros(length(sample_sizes))

for (i, sample_size) in enumerate(sample_sizes)
    println(i, " out of ", length(sample_sizes))

    data = DeltaData(
        params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
            :x2 => Uniform(-3.14159265359, 3.14159265359),
            :x3 => Uniform(-3.14159265359, 3.14159265359),
            :x4 => Uniform(-3.14159265359, 3.14159265359)),
        N = sample_size
    )

    X = GlobalSensitivityAnalysis.sample(data)
    Y = ishigami(X)
    results = analyze(data, X, Y; num_resamples = 500)
    
    fig3_biased_deltas[i] = results[:delta_biased][4]
    fig3_deltas[i] = results[:delta][4]
    fig3_deltas_CI[i] = results[:delta_conf][4]
end

df = DataFrame(
    samples = sample_sizes,
    biased_deltas = fig3_biased_deltas, 
    deltas = fig3_deltas,
    CI_high = fig3_deltas + fig3_deltas_CI,
    CI_low = fig3_deltas - fig3_deltas_CI
)

p = df |> @vlplot(
        x = {field = :samples, title = "sample size"},
        y = {title = "delta sensitivity measure"},
        width = 500, 
        height = 250
    ) + 
    @vlplot(
        mark={
            type = :line, 
            color = :darkorange,
            point={filled = false, fill = :white, color = :darkorange}
        },
        y=:biased_deltas
    ) + 
    @vlplot(
        mark={
            :line, 
            point={filled = false, fill = :white}
        },
        y=:deltas
    ) +
    @vlplot(
        mark=:area,
        y = :CI_high,
        y2 = :CI_low,
        opacity={value=0.25}
    )

save(joinpath(output_dir, "ReplicateFig3.png"), p)
