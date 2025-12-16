# DEQuadrature.jl

[![CI](https://github.com/JuliaReliab/DEQuadrature.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaReliab/DEQuadrature.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaReliab/DEQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReliab/DEQuadrature.jl)

Fast, robust numerical integration on finite and semi-infinite intervals using the Double Exponential (tanh–sinh) quadrature.

DEQuadrature.jl provides a simple, allocation‑light API to integrate 1‑D functions with good accuracy for smooth and rapidly decaying integrands.

## Installation

Once registered in the General registry:

```julia
using Pkg
Pkg.add("DEQuadrature")
```

Until then, you can install from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaReliab/DEQuadrature.jl.git")
```

## Quick start

```julia
using DEQuadrature

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
      maxiter::Int = 12)
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

### Arbitrary precision (BigFloat)

DEQuadrature also supports arbitrary precision via a generic method:

```julia
deint(f, lower::T, upper::T;
        reltol::T = T(1.0e-9),
        abstol::T = eps(T),
        d::Int = 8,
        maxiter::Int = 12) where {T<:Real}
```

Usage with `BigFloat`:

```julia
using DEQuadrature

setprecision(256) do
      # ∫_0^∞ e^{-x} dx = 1
      res = deint(x -> exp(-x), BigFloat(0), BigFloat(Inf); reltol=BigFloat(1e-20))
      @show res.s  # ≈ 1 with high precision

      # Finite interval [a,b]
      a = BigFloat(-2); b = BigFloat(3)
      res2 = deint(x -> one(BigFloat), a, b; reltol=BigFloat(1e-18))
      @show res2.s  # ≈ b - a
end
```

Tips:
- Set `setprecision` to suit your needs and tighten `reltol/abstol` accordingly.
- For very tight tolerances or slowly decaying integrands, consider increasing `maxiter` and/or the initial divisions `d`.

### Tuning for performance

The default parameters `(d=8, maxiter=12)` balance speed and accuracy for most smooth, well-behaved integrands. Adjust based on your needs:

**maxiter** (default: 12)
- Controls the refinement depth. Each iteration doubles the node count.
- Sufficient for `reltol ~1e-8` with `Float64`. For tighter tolerances, increase modestly (e.g., to 14–16).
- Be cautious: `maxiter=16` can require ~262k nodes; prefer tightening `reltol`/`abstol` first.

**d** (default: 8)
- Initial node count per interval. Larger `d` front-loads computation but may need fewer refinement steps.
- Try `d=10–16` for smooth integrands; `d=6–8` for rougher ones or memory-constrained settings.

**reltol/abstol** (defaults: `1e-8`/`eps(Float64)`)
- Convergence is based on the relative/absolute error between iterations.
- Tighten these *first* for higher accuracy; increasing `maxiter` is a fallback when tolerance targets are nearly met but time/memory is less critical.

Advanced users can also construct the DE mapping for an arbitrary `T` via:

```julia
deformula_zero_to_inf(T)       # semi-infinite (0, ∞)
deformula_minus_one_to_one(T)  # finite interval (-1, 1)
```
and call the internal `_deint(f, formula::Formula{T}; ...)` directly if needed.

### Exported symbols

- `deint` — high‑level integration API
- `deformulaZeroToInf`, `deformulaMinusOneToOne` — prebuilt DE mappings (advanced use)
      - For arbitrary precision: use `deformula_zero_to_inf(T)` / `deformula_minus_one_to_one(T)`

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
