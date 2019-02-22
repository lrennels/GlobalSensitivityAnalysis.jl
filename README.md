[![Build Status](https://travis-ci.org/lrennels/GlobalSensitivityAnalysis.jl.svg?branch=master)](https://travis-ci.org/lrennels/GlobalSensitivityAnalysis.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/github/lrennels/GlobalSensitivityAnalysis.jl?branch=master&?svg=true)](https://ci.appveyor.com/project/lrennels/GlobalSensitivityAnalysis.jl/branch/master)
[![codecov](https://codecov.io/gh/lrennels/GlobalSensitivityAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lrennels/GlobalSensitivityAnalysis.jl)
[![Coverage Status](https://img.shields.io/coveralls/github/lrennels/GlobalSensitivityAnalysis.jl/master.svg)](https://coveralls.io/github/lrennels/GlobalSensitivityAnalysis.jl?branch=master)

# Global Sensitivity Analysis

A Julia package originally based on [SALib](https://github.com/SALib/SALib) (Herman and Usher, 2017) which implements global sensitivity analysis methods.

The package currently includes the following methods:

- Sobol Sensitivity Analysis ([Sobol 2001](http://www.sciencedirect.com/science/article/pii/S0378475400002706), [Saltelli 2002](http://www.sciencedirect.com/science/article/pii/S0010465502002801), [Saltelli et al. 2010](http://www.sciencedirect.com/science/article/pii/S0010465509003087))

**References**

```
Herman, J. and Usher, W. (2017) SALib: An open-source Python library for sensitivity analysis. 
Journal of Open Source Software, 2(9).

Saltelli, A. (2002).  "Making best use of model evaluations to compute sensitivity indices." 
Computer Physics Communications,145(2):280-297, doi:10.1016/S0010-4655(02)00280-1.

Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and S. Tarantola (2010).  
"Variance based sensitivity analysis of model output.  Design and estimator for the total 
sensitivity index." Computer Physics Communications, 181(2):259-270, 
doi:10.1016/j.cpc.2009.09.018.

Sobol, I. M. (2001).  "Global sensitivity indices for nonlinear mathematical models and their 
Monte Carlo estimates."  Mathematics and Computers in Simulation, 55(1-3):271-280, 
doi:10.1016/S0378-4754(00)00270-6.
```

## Copyright Info

Some of the code in this package is derivative code from:

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

The code in `sobol_sequence.jl` is derivative code of `SALib` as mentioned above, 
which in turn derived the code from a source with the following copyright:

    Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
