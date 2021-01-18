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

# load SALib values
SALib_deltas = [0.03675842, 0.03512737, 0.04431219, 0.03533525, 0.0386379 , 0.02777811]
SALib_deltas_CI = [0.02164714, 0.01424834, 0.01058716, 0.0073074 , 0.00521382, 0.00358639]

# get my values
sample_sizes = 2 .^ [i for i in 9:14]
GSA_biased_deltas = zeros(length(sample_sizes))
GSA_deltas = zeros(length(sample_sizes))
GSA_deltas_CI = zeros(length(sample_sizes))

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
    results = analyze(data, X, Y; num_resamples = 1000)
    
    GSA_biased_deltas[i] = results[:delta_biased][4]
    GSA_deltas[i] = results[:delta][4]
    GSA_deltas_CI[i] = results[:delta_conf][4]
end

# build dataframe
df = DataFrame(
    samples = sample_sizes,
    GSA_biased_deltas = GSA_biased_deltas, 
    GSA_deltas = GSA_deltas,
    GSA_CI_high = GSA_deltas + GSA_deltas_CI,
    GSA_CI_low = GSA_deltas - GSA_deltas_CI,
    SALib_deltas = SALib_deltas,
    SALib_CI_high = SALib_deltas + SALib_deltas_CI,
    SALib_CI_low = SALib_deltas - SALib_deltas_CI,
)

# plot
p = df |> @vlplot(
        x = {field = :samples, title = "sample size"},
        width = 500, 
        height = 250
    ) + 
    @vlplot(
        mark={
            type = :line, 
            color = :darkorange,
            point={filled = false, fill = :white, color = :darkorange}
        },
        y=:GSA_biased_deltas
    ) + 
    @vlplot(
        mark=:area,
        y = :GSA_CI_high,
        y2 = :GSA_CI_low,
        opacity={value=0.45}
    ) +
    @vlplot(
        mark= {:area, color = :darkgreen},
        y = {field = :SALib_CI_high, title = "delta sensitivity measure"},
        y2 = {field = :SALib_CI_low},
        opacity={color = :darkgreen, value=0.15}
    ) +
    @vlplot(
        mark={
            :line, 
            color = :darkgreen,
            point={filled = false, fill = :white, color = :darkgreen}
        },
        y=:SALib_deltas
    ) +
    @vlplot(
        mark={
            :line, 
            point={filled = false, fill = :white}
        },
        y=:GSA_deltas
    ) 
    

save(joinpath(output_dir, "Replicate_Fig3_resample1000.png"), p)
