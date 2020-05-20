
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

function _calcWeight!(data::Vector{Tuple{T,T,T}}, t::T, f, phi, phidash; abstol::T = eps(T)) where {T <: Real}
    local xtmp::T = phi(t)
    local wtmp::T = phidash(t) * f(xtmp)
    if !isnan(wtmp) && wtmp > abstol
        if !isfinite(wtmp)
            error("Error: weight becomes NaN")
        end
        push!(data, (t, xtmp, wtmp))
    end
end

"""
deint(f, formula; reltol::T = 1.0e-8, abstol::T = eps(T), d = 8, maxiter = 31)

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
- the value of integration
- a sequence for divides on the transformed domain
- a sequence for divides; t_i
- a sequence of weights (unscaled); w_i
- a scale parameter for w_i; h
"""

function deint(f, formula::Formula{T};
        reltol::T = 1.0e-8, abstol::T = eps(T), d = 8, maxiter = 31) where {T <: Real}
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
        prev = s
        if iter > maxiter
            error("Error: The number of iteration attains maxtier")
        end
        h /= 2
        for t = LinRange(lower+h, upper-h, d)
            _calcWeight!(data, t, f, formula.phi, formula.phidash, abstol=abstol)
        end
        d *= 2
        s = sum([x[3] for x in data]) * h
        aerror = s - prev
        rerror = aerror / prev
        if abs(rerror) < reltol
            break
        end
    end
    sort!(data, by=first)
    (s, map(x->x[1], data), map(x->x[2], data), map(x->x[3], data), h)
end
