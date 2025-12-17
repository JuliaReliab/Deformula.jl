
export deformulaMinusOneToOne, deformulaZeroToInf, deint

struct Formula{T <: Real}
    range::Tuple{T,T}
    phi
    phidash
end

# (Float64 constant kept via constructor below for backward compatibility)
"""
    deformula_zero_to_inf(::Type{T}=Float64) where {T<:Real}

Constructor for the semi-infinite range DE formula specialized to numeric type `T`.
"""
deformula_zero_to_inf(::Type{T}=Float64) where {T<:Real} = Formula{T}( 
    (-T(6.8), T(6.8)),
    t -> exp(T(pi) * sinh(t) / T(2)),
    t -> T(pi) * cosh(t) * exp(T(pi) * sinh(t) / T(2)) / T(2)
)

const deformulaZeroToInf = deformula_zero_to_inf(Float64)

# (Float64 constant kept via constructor below for backward compatibility)
"""
    deformula_minus_one_to_one(::Type{T}=Float64) where {T<:Real}

Constructor for the finite range [-1, 1] DE formula specialized to numeric type `T`.
"""
deformula_minus_one_to_one(::Type{T}=Float64) where {T<:Real} = Formula{T}(
    (-T(3.0), T(3.0)),
    t -> tanh(T(pi) * sinh(t) / T(2)),
    t -> begin
        pisinh2 = T(pi) * sinh(t) / T(2)
        sech = one(T) / cosh(pisinh2)
        T(pi) * cosh(t) * sech * sech / T(2)
    end
)

const deformulaMinusOneToOne = deformula_minus_one_to_one(Float64)

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
    _deint(f, formula; reltol::T = 1.0e-8, abstol::T = eps(T), d = 8, maxiter = 12)

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
        reltol::T = 1.0e-9, abstol::T = eps(T), d = 8, maxiter = 12) where {T <: Real}
    local lower::T = formula.range[1]
    local upper::T = formula.range[2]
    local h::T = (upper - lower) / d

    # Pre-allocate with an upper-bound estimate to reduce reallocations
    data = Vector{Tuple{T,T,T}}()
    if maxiter >= 1
        # Total appended points ≈ (d+1) + d * (1 + 2 + ... + 2^(maxiter-2)) = (d+1) + d * (2^(maxiter-1) - 1)
        local expected_len::Int = d + 1 + d * ((1 << (maxiter - 1)) - 1)
        sizehint!(data, expected_len)
    else
        sizehint!(data, d + 1)
    end
    
    for t = LinRange(lower, upper, d+1)
        _calcWeight!(data, t, f, formula.phi, formula.phidash, abstol=abstol)
    end
    
    # Compute sum directly without intermediate array
    local s::T = zero(T)
    for item in data
        s += item[3]
    end
    s *= h
    
    local iter = 1
    local prev::T = s
    while true
        iter += 1
        if iter > maxiter
            @warn "Max iterations reached in _deint; returning last estimate" maxiter=maxiter s=s h=h d=d
            break
        end
        h /= 2
        for t = LinRange(lower+h, upper-h, d)
            _calcWeight!(data, t, f, formula.phi, formula.phidash, abstol=abstol)
        end
        d *= 2
        
        # Compute sum directly
        s = zero(T)
        for item in data
            s += item[3]
        end
        s *= h
        
        aerror = abs(s - prev)
        denom = max(abs(prev), one(T))
        rerror = aerror / denom
        if (aerror <= abstol) || (rerror <= reltol)
            break
        end
        prev = s
    end
    sort!(data, by=first)
    
    # Pre-allocate output vectors
    n = length(data)
    t_vec = Vector{T}(undef, n)
    x_vec = Vector{T}(undef, n)
    w_vec = Vector{T}(undef, n)
    @inbounds for i in 1:n
        t_vec[i] = data[i][1]
        x_vec[i] = data[i][2]
        w_vec[i] = data[i][3]
    end
    
    (s=s, t=t_vec, x=x_vec, w=w_vec, h=h)
end

"""
    deint(f, lower, upper; ...)

Apply the double exponential (tanh-sinh) procedure to `f` on [`lower`, `upper`].

This function returns not only the integral value, but also the nodes and
weights used in the computation. It is intended for research and analysis
purposes rather than as a black-box integrator.

Returns a named tuple `(s, t, x, w, h)` with the integral estimate, nodes in the
transformed domain, mapped points, weights, and scale factor:
- `s`: integral estimate
- `t`: nodes in the transformed domain
- `x`: mapped nodes in the original domain
- `w`: unscaled weights corresponding to `x`
- `h`: scale factor applied to weights
"""
function deint(f, lower::T, upper::T;
    reltol::T=T(1.0e-9), abstol::T=eps(T), d=8, maxiter=12) where {T<:Real}

    if isinf(upper) && upper > zero(T)
        if lower == zero(T)
            return _deint(f, deformula_zero_to_inf(T), reltol=reltol, abstol=abstol, d=d, maxiter=maxiter)
        else
            result = _deint(deformula_zero_to_inf(T), reltol=reltol, abstol=abstol, d=d, maxiter=maxiter) do x
                f(x + lower)
            end
            x = result.x .+ lower
            return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
        end
    else
        if lower == -one(T) && upper == one(T)
            return _deint(f, deformula_minus_one_to_one(T), reltol=reltol, abstol=abstol, d=d, maxiter=maxiter)
        else
            d_half = (upper - lower) / T(2)
            result = _deint(deformula_minus_one_to_one(T), reltol=reltol, abstol=abstol, d=d, maxiter=maxiter) do x
                f(d_half * (x + one(T)) + lower) * d_half
            end
            x = @. d_half * (result.x + one(T)) + lower
            return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
        end
    end

end
"""
        deint(f, lower, upper; ...)

Apply the double exponential (tanh-sinh) procedure to `f` on [`lower`, `upper`].

This function returns not only the integral value, but also the nodes and
weights used in the computation. It is intended for research and analysis
purposes rather than as a black-box integrator.

Returns a named tuple `(s, t, x, w, h)` with the integral estimate, nodes in the
transformed domain, mapped points, weights, and scale factor:
- `s`: integral estimate
- `t`: nodes in the transformed domain
- `x`: mapped nodes in the original domain
- `w`: unscaled weights corresponding to `x`
- `h`: scale factor applied to weights
"""
function deint(f, lower::Float64, upper::Float64;
    reltol::Float64=1.0e-8, abstol::Float64=eps(Float64), d=8, maxiter=12)::NamedTuple{(:s, :t, :x, :w, :h), Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64}}
    
    if isinf(upper) && upper > 0
        # Handle [0, Inf) or [lower, Inf) cases
        if lower == 0.0
            return _deint(f, deformulaZeroToInf, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter)
        else
            result = _deint(deformulaZeroToInf, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter) do x
                f(x + lower)
            end
            # Use broadcasting for efficiency
            x = result.x .+ lower
            return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
        end
    else
        # Handle finite interval cases
        if lower == -1.0 && upper == 1.0
            return _deint(f, deformulaMinusOneToOne, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter)
        else
            # General finite interval
            d_half = (upper - lower) / 2
            result = _deint(deformulaMinusOneToOne, reltol=reltol, abstol=abstol, d=d, maxiter=maxiter) do x
                f(d_half * (x + 1) + lower) * d_half
            end
            # Use broadcasting and fused operations for efficiency
            x = @. d_half * (result.x + 1.0) + lower
            return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
        end
    end
end
