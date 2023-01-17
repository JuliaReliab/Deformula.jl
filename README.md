# Deformula

[![Build Status](https://travis-ci.com/okamumu/Deformula.jl.svg?branch=master)](https://travis-ci.com/okamumu/Deformula.jl)
[![Codecov](https://codecov.io/gh/okamumu/Deformula.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/okamumu/Deformula.jl)
[![Coveralls](https://coveralls.io/repos/github/okamumu/Deformula.jl/badge.svg?branch=master)](https://coveralls.io/github/okamumu/Deformula.jl?branch=master)

Deformula.jl provides numerical quadrature of functions of one variable over a finite or infinite interval with double exponential formulas.

## Installation

This is not in the official Julia package yet. Please run the following command to install it.
```
using Pkg; Pkg.add(PackageSpec(url="https://github.com/JuliaReliab/Deformula.jl.git"))
```

## Load module

Load the module:
```
using Deformula
```

## How to use

The package provides a function `deint` to obtain numerical quadrature.

```julia
deint(f, lower, upper; reltol::T = 1.0e-8, abstol::T = eps(T), d = 8, maxiter = 16)
```

This function computes the numerical integration for `f` on `[lower, upper]` with double exponential formula.

$$
\int f(x) dx \approx h * \sum_i w_i
$$

- Parameters:
    - f: integrand function
    - lower: lower value of range. This should be finite.
    - upper: upper value of range. This is allowed to take infinite.
    - reltol: tolerance for relative errors
    - abstol: tolerance for absolute errors
    - d: the initial number of divides. The default is `d=8`
    - maxiter: the maximum number of iterations to increase the number of divides twice.
- Return value:
    - s: the value of integration
    - t: a sequence for divides on the transformed domain
    - x: a sequence for divides; t_i
    - w: a sequence of weights (unscaled); w_i
    - h: a scale parameter for w_i; h

## Example

Obtain the numerical quadrature of $f(x) = e^{-x}$ on $[0, \infty)$. 

```julia
result = deint(x -> exp(-x), 0.0, Inf) ## lower and upper should be Float64.
println(result.s)
println(result.h * sum(result.w))
```

## Note

We use the integartion with semi-infinite interval and finite interval.
