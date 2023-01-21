
export deformulaMinusOneToOne, deformulaZeroToInf, deint

struct Formula{T <: Real}
    range::Tuple{T,T}
    phi
    phidash
end

const deformulaZeroToInf = Formula(
    (-6.8, 6.8),
    t -> exp(pi * sinh(t) / 2),
    t -> pi * cosh(t) * exp(pi * sinh(t) / 2) / 2
)

const deformulaMinusOneToOne = Formula(
    (-3.0, 3.0),
    t -> tanh(pi * sinh(t) / 2),
    t -> pi * cosh(t) * (1 / cosh(pi * sinh(t) / 2)) * (1 / cosh(pi * sinh(t) / 2)) / 2
)

# struct DeformulaResult{T <: Real}
#     s::T
#     t::Vector{T}
#     x::Vector{T}
#     w::Vector{T}
#     h::T
# end

function _calcWeight!(data::Vector{Tuple{T,T,T}}, t::T, f, phi, phidash; abstol::T = eps(T)) where {T <: Real}
    local xtmp::T = phi(t)
    local wtmp::T = phidash(t) * f(xtmp)
    if !isnan(wtmp) && abs(wtmp) > abstol
        if !isfinite(wtmp)
            error("Error: weight becomes NaN")
        end
        push!(data, (t, xtmp, wtmp))
    end
end

"""
_deint(f, formula; reltol::T = 1.0e-8, abstol::T = eps(T), d = 8, maxiter = 16)

Compute the numerical integration for f with double exponential formula.

    int f(x) dx = h * sum w_i

Parameters:
- f: integrand function
- formula: double expoential formula. 
- reltol: tolerance for relative errors
- abstol: tolerance for absolute errors
- d: the initial number of divides
- maxiter: the maximum number of iterations to increase the number of divides twice.
Return value (tuple):
- s: the value of integration
- t: a sequence for divides on the transformed domain
- x: a sequence for divides; x_i
- w: a sequence of weights (unscaled); w_i
- h: a scale parameter for w_i
"""

function _deint(f, formula::Formula{T};
        reltol::T = 1.0e-9, abstol::T = eps(T), d = 8, maxiter = 16) where {T <: Real}
    local lower::T = formula.range[1]
    local upper::T = formula.range[2]
    local h::T = (upper - lower) / d

    data = Vector{Tuple{T,T,T}}()
    for t = LinRange(lower, upper, d+1)
        _calcWeight!(data, t, f, formula.phi, formula.phidash, abstol=abstol)
    end
    local s::T = sum([x[3] for x in data]) * h
    
    local iter = 1
    while true
        iter += 1
        prev = s + 1.0
        if iter > maxiter
            error("Warning: The number of iteration attains maxtier $(maxiter)")
        end
        h /= 2
        for t = LinRange(lower+h, upper-h, d)
            _calcWeight!(data, t, f, formula.phi, formula.phidash, abstol=abstol)
        end
        d *= 2
        s = sum([x[3] for x in data]) * h
        aerror = s + 1.0 - prev
        rerror = aerror / prev
        if abs(rerror) < reltol
            break
        end
    end
    sort!(data, by=first)
    (s=s, t=map(x->x[1], data), x=map(x->x[2], data), w=map(x->x[3], data), h=h)
end

"""
    deint(f, lower::Float64, upper::Float64;
        reltol::Float64=1.0e-8, abstol::Float64=eps(Float64), d=8, maxiter=16)

Compute the numerical integration for `f` on [`lower`, `upper`] with double exponential formula.

    int f(x) dx = h * sum w_i

# Arguments
  - `f`: integrand function
  - `lower::Float64`: lower value of range. This should be finite.
  - `upper::Float64`: upper value of range. This is allowed to take infinite.
  - `reltol::Float64=1.0e-8`: tolerance for relative errors
  - `abstol::Float64=eps(Float64)`: tolerance for absolute errors
  - `d=8`: the initial number of divides
  - `maxiter=16`: the maximum number of iterations to increase the number of divides twice.

# Return values
  - `s`: the value of integration
  - `t`: a sequence for divides on the transformed domain
  - `x`: a sequence for divides; t_i
  - `w`: a sequence of weights (unscaled); w_i
  - `h`: a scale parameter for w_i; h

# Examples
```julia-repl
julia> f(x) = exp(x);
julia> result = deint(f, 0.0, 1.0);
```

```julia-repl
julia> result = deint(0.0, Inf64) do x
    2.0*exp(-2.0*x)
end
```
"""
function deint(f, lower::Float64, upper::Float64;
    reltol::Float64=1.0e-8, abstol::Float64=eps(Float64), d=8, maxiter=16)::NamedTuple{(:s, :t, :x, :w, :h), Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64}}
    if lower == 0.0 && upper == Inf64
        _deint(f, deformulaZeroToInf, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter)
    elseif lower != 0.0
        result = _deint(deformulaZeroToInf, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter) do x
            f(x + lower)
        end
        x = result.x .+ lower
        (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
    elseif lower == -1.0 && upper == 1.0
        _deint(f, deformulaMinusOneToOne, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter)
    elseif upper != Inf64
        result = _deint(deformulaMinusOneToOne, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter) do x
            d = (upper - lower)/2
            f(d*(x + 1) + lower) * d
        end
        x = ((upper - lower)/2) .* (result.x .+ 1.0) .+ lower
        (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
    else
        error("lower and upper should be finite or semi-infinite.")
    end
end
