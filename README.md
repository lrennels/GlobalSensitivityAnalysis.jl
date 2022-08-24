![](https://github.com/lrennels/GlobalSensitivityAnalysis.jl/workflows/Run%20CI%20on%20master/badge.svg)
[![codecov](https://codecov.io/gh/lrennels/GlobalSensitivityAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lrennels/GlobalSensitivityAnalysis.jl)
[![Coverage Status](https://coveralls.io/repos/github/lrennels/GlobalSensitivityAnalysis.jl/badge.svg?branch=master)](https://coveralls.io/github/lrennels/GlobalSensitivityAnalysis.jl?branch=master)

# Global Sensitivity Analysis

A Julia package which implements global sensitivity analysis methods.

The package currently includes the following methods:

- Sobol Sensitivity Analysis ([Sobol 2001](http://www.sciencedirect.com/science/article/pii/S0378475400002706), [Saltelli 2002](http://www.sciencedirect.com/science/article/pii/S0010465502002801), [Saltelli et al. 2010](http://www.sciencedirect.com/science/article/pii/S0010465509003087))

- Delta Moment-Independent Measure ([Borgonovo 2007](http://www.sciencedirect.com/science/article/pii/S0951832006000883), [Plischke et al. 2013](http://www.sciencedirect.com/science/article/pii/S0377221712008995))

## SALib

Much of this package is based on/transalted from[SALib](https://github.com/SALib/SALib) (Herman and Usher, 2017) which implements several global sensitivity analysis measures in Python.  The present package seeks to implement several of these same algorithms in Julia along with providing a clear, user-friendly API.

The license terms are provided below:

The MIT License (MIT)
Copyright (c) 2013-2017 Jon Herman, Will Usher, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## The API

The API contains two primary functions: `sample` and `analyze`. These two functions call methods based on the type parameterization of their `data` argument, which is either of type `SobolData` or `DeltaData`.  

**Note:** For now the `sample` function will call the most used sampling protocol for the particular method, (Sobol sequence for Sobol method and Latin Hypercube sampling for Delta method), but in this future this will be rearranged and generalized since, for example, the Delta method can also just as well use Sobol sequence sampling and other methods.

### Sobol Sensitivity Analyis

Sampling with `sample` is the first of the two main steps in an analysis, generating the model inputs to be run through a model of choice and produce the outputs analyzed in the `analyze` function.  The signature for this function is as follows.

```julia
    sample(data::SobolData)

Generate a matrix containing the model inputs for Sobol sensitivity analysis with 
the information in the `data`. In this function we apply Saltelli's 
extension of the Sobol  sequence. Saltelli's scheme extends the Sobol sequence in 
a way to reduce the error rates in the resulting sensitivity index calculations. 
The resulting matrix has `N` * (`D` + 2) rows, where `D` is the number of parameters 
and `N` is the number of samples.
```

The single argument to this function is of type `SobolData`, a custom type designed to hold all information needed for sampling and analysis. A `SobolData` struct is parameterized by a `params` dictionary (NOTE that this must be an `OrderedDict`, not a `Dict`) which maps parameter names to their Distributions, `calc_second_order` determining whether or not to calculate second-order sensitivity indices, and the desired number of runs `N`.

```julia
    SobolData

A struct which holds all information needed for the sampling and analysis of a
specific problem using Sobol Analysis:

`params::Union{OrderedDict{Symbol, <:Any}, Nothing} = nothing`: a dictionary mapping parameter names to their Distribution
`calc_second_order::Bool = true`: whether or not to calculate second order sensitivity indices
`N::Int = 1000`: the number of runs
```

After sampling with `sample`, use the resulting matrix of parameter combinations to run your model, producing a vector of results.  The next and final step is to analyze the results with your `model_output` using the `analyze` function with the signature below. This function takes the same `SobolData` as `sample`, as well as the `model_output` vector and produces a dictionary of results.  This dictionary will include the `:firstorder`, `:totalorder` indices and (optionally) confidence intervals for each parameter.

```julia
    function analyze(data::SobolData, model_output::AbstractArray{<:Number, S}; num_resamples::Union{Nothing, Int} = 1_000, conf_level::Union{Nothing, Number} = 0.95, progress_meter::Bool = true, N_override::Union{Nothing, Integer}=nothing) 

Performs a Sobol Analysis on the `model_output` produced with the problem defined by the information in `data` and returns the a dictionary of results with the sensitivity indices and respective confidence intervals for each of the
parameters defined using the `num_resamples` and `conf_level` keyword args. If these are Nothing than no confidence intervals will be calculated. The `progress_meter` keyword argument indicates whether a progress meter will be displayed and defaults to true. The `N_override` keyword argument allows users to override the `N` used in a specific `analyze` call to analyze just a subset (useful for convergence graphs).
```

An example of the basic flow can be found in `src/main.jl` using the Ishigami test function in `src/test_functions/ishigami.jl`, and is copied and commented below for convenience.

```julia
using Distributions
using DataStructures

include("sample_sobol.jl")
include("analyze_sobol.jl")
include("test_functions/ishigami.jl")

# define the data
data = SobolData(
    params = OrderedDict(:x1 => Uniform(-3.14159265359, 3.14159265359),
        :x2 => Uniform(-3.14159265359, 3.14159265359),
        :x3 => Uniform(-3.14159265359, 3.14159265359)),
    N = 1000
)

# generate samples using Sobol sequence
samples = sample(data)

# run model (example)
Y = ishigami(samples)

# perform Sobol Analysis
analyze(data, Y)
```
### Delta Moment-Independent Measure

Sampling with `sample` is the first of the two main steps in an analysis, generating the model inputs to be run through a model of choice and produce the outputs analyzed in the `analyze` function.  The signature for this function is as follows.

```julia
    sample(data::DeltaData)

Generate a matrix containing the model inputs for Delta Moment-Independent Measure sensitivity analysis with the information in the `data`. In this function we apply Latin Hypercube Sampling. The resulting matrix has `N` columns * `D` rows, where `D` is the number of parameters and `N` is the number of samples.
```

The single argument to this function is of type `DeltaData`, a custom type designed to hold all information needed for sampling and analysis. A `DeltaData` struct is parameterized by a `params` dictionary (NOTE that this must be an `OrderedDict`, not a `Dict`) which maps parameter names to their Distributions and the desired number of runs `N`.

```julia
    DeltaData

A struct which holds all information needed for the sampling and analysis of a
specific problem using Sobol Analysis:

`params::Union{OrderedDict{Symbol, <:Any}, Nothing} = nothing`: a dictionary mapping parameter names to their Distribution
`N::Int = 1000`: the number of runs
```

After sampling with `sample`, use the resulting matrix of parameter combinations to run your model, producing a vector of results.  The next and final step is to analyze the results with your `model_output` using the `analyze` function with the signature below. This function takes the same `DeltaData` as `sample`, as well as the `model_output` vector and produces a dictionary of results.  This dictionary will include the `:firstorder`, `:delta` indices and confidence intervals for each parameter.

```julia
    function analyze(data::DeltaData, model_input::AbstractArray{<:Number, S1}, model_output::AbstractArray{<:Number, S2}; num_resamples::Int = 1_000, conf_level::Number = 0.95, progress_meter::Bool = true, N_override::Union{Nothing, Integer}=nothing)

Performs a Delta Moment-Independent Analysis on the `model_output` produced with  the problem defined by the information in `data` and `model_input` and returns a dictionary of results with the sensitivity indices and respective confidence  intervals for each of the parameters defined using the `num_resamples` and  `conf_level` keyword args.  The `progress_meter` keyword argument indicates whether a  progress meter will be displayed and defaults to true. The `N_override` keyword  argument allows users to override the `N` used in a specific `analyze` call to analyze just a subset (useful for convergence graphs).
```

## References

References from the peer-reviewed literature include:

    Borgonovo, E. (2007). A new uncertainty importance measure. Reliability Engineering 
    & System Safety, 92(6), 771-784. Chicago.

    Herman, J. and Usher, W. (2017) SALib: An open-source Python library for sensitivity 
    analysis. Journal of Open Source Software, 2(9).

    Plischke, E., Borgonovo, E., & Smith, C. L. (2013). Global sensitivity measures 
    from given data. European Journal of Operational Research, 226(3), 536-550.

    Saltelli, A. (2002).  "Making best use of model evaluations to compute sensitivity indices." 
    Computer Physics Communications,145(2):280-297, doi:10.1016/S0010-4655(02)00280-1.

    Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and S. Tarantola (2010).  
    "Variance based sensitivity analysis of model output.  Design and estimator 
    for the total sensitivity index." Computer Physics Communications, 181(2):259-270, 
    doi:10.1016/j.cpc.2009.09.018.

    Sobol, I. M. (2001).  "Global sensitivity indices for nonlinear mathematical 
    models and their Monte Carlo estimates."  Mathematics and Computers in Simulation, 
    55(1-3):271-280, doi:10.1016/S0378-4754(00)00270-6.

## Copyright Information

Some of the code in this package is derivative code based on the work of John Herman, Will Usher, and others:

    The MIT License (MIT)

    Copyright (c) 2013-2017 Jon Herman, Will Usher, and others.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
