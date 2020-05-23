# QuadGG.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/tbeason/QuadGG.jl.svg?branch=master)](https://travis-ci.com/tbeason/QuadGG.jl)
[![codecov.io](http://codecov.io/github/tbeason/QuadGG.jl/coverage.svg?branch=master)](http://codecov.io/github/tbeason/QuadGG.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://tbeason.github.io/QuadGG.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://tbeason.github.io/QuadGG.jl/dev)
-->

`QuadGG.jl` provides a pure Julia implementation of the MATLAB quadrature routines from Gander & Gauschi (2000) ["Adaptive Quadrature -- Revisited"](ftp://ftp.inf.ethz.ch/pub/publications/tech-reports/3xx/306.ps.gz). 

The first is an adaptive Gauss-Lobatto method, while the second is an adaptive Simpson's rule. The Gauss-Lobatto routine seems to be the preferred method of the two.

## Usage

There are two exported functions: `adaptlob` for the adaptive Gauss-Lobatto routine and `adaptsim` for the adaptive Simpson's rule.

```julia
using QuadGG

adaptlob(x->x^2,0,1) #  == 1/3
adaptlob(sqrt,0,1;tol=1e-16) #  == 2/3

adaptsim(x->x^2,0,1) #  == 1/3
adaptsim(sqrt,0,1;tol=1e-16) #  == 2/3
```

## Alternative Packages

I wrote this package (ported the code is more like it) so that I could take advantage of a particular feature of the problem I was solving. [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) and [HCubature.jl](https://github.com/JuliaMath/HCubature.jl) offer additional features compared to this library and are likely better choices for general numerical integration problems. 

The biggest conceptual difference from `QuadGK.jl` is that the Gauss-Lobatto scheme does evaluate the integral at the interval endpoints and they use different adaptive sampling techniques. The `QuadGK.jl` methods are very lean, therefore if your integrand is _really_ smooth they can integrate pretty much instantly. This package would use many more function evaluations in those cases. For more difficult integrals, the speed differences between the package are not so great (and sometimes this package is faster).
