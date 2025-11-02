# Deformula.jl

[![CI](https://github.com/JuliaReliab/Deformula.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaReliab/Deformula.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaReliab/Deformula.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReliab/Deformula.jl)

Fast, robust numerical integration on finite and semi-infinite intervals using the Double Exponential (tanh–sinh) quadrature.

Deformula.jl provides a simple, allocation‑light API to integrate 1‑D functions with good accuracy for smooth and rapidly decaying integrands.

## Installation

Once registered in the General registry:

```julia
using Pkg
Pkg.add("Deformula")
```

Until then, you can install from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaReliab/Deformula.jl.git")
```

## Quick start

```julia
using Deformula

# Integrate e^{-x} on [0, ∞)
res = deint(x -> exp(-x), 0.0, Inf64)
@show res.s           # integral value
@show res.h * sum(res.w)  # scaled weights (sanity check)

# Integrate on a finite interval [a,b]
a, b = -2.5, 3.7
res2 = deint(x -> 1.0, a, b)
@show res2.s  # ≈ b - a
```

## API

```julia
deint(f, lower::Float64, upper::Float64;
      reltol::Float64 = 1.0e-8,
      abstol::Float64 = eps(Float64),
      d::Int = 8,
      maxiter::Int = 16)
```

Computes

$$
\int_{lower}^{upper} f(x)\,dx \;\approx\; h\,\sum_i w_i
$$

Return value is a named tuple `(s, t, x, w, h)`:

- `s`: integral estimate
- `t`: nodes in the transformed domain
- `x`: nodes in the original domain
- `w`: unscaled weights
- `h`: scale for weights (so that `h*sum(w) ≈ s`)

Notes:

- Convergence uses combined absolute/relative checks. Increase `maxiter` or `d` if you need tighter accuracy.
- On reaching `maxiter` before convergence, a warning is emitted and the last estimate is returned.

### Exported symbols

- `deint` — high‑level integration API
- `deformulaZeroToInf`, `deformulaMinusOneToOne` — prebuilt DE mappings (advanced use)

## Compatibility

- Julia 1.x (see `Project.toml` for details)

## Development

- Run tests:

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

Contributions are welcome. Please open an issue or pull request.

## License

This project is licensed under the terms of the LICENSE file included in this repository.
