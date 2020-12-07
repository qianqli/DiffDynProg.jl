#=
Created on 07/12/2020 09:23:01
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the smooth max operators and their gradients.
=#

using LinearAlgebra: ⋅, norm
using StatsBase: mean
using ChainRulesCore
import ChainRulesCore: frule, rrule
import Base: min, max, minimum, maximum

# Type hierarchy
# --------------

abstract type MaxOperator end

# ordinary max
struct Max <: MaxOperator end

# leaks a fraction p
struct LeakyMax{T<:AbstractFloat} <: MaxOperator
    p::Float64
    function LeakyMax(p::T=0.1) where {T<:AbstractFloat}
        @assert 0.0 < p < 1.0 "`p` has to be in ]0, 1["
        return new{T}(p)
    end
end

struct EntropyMax{T<:AbstractFloat} <: MaxOperator
    γ::T
    function EntropyMax(γ::T=1.0) where {T<:AbstractFloat}
        @assert γ > 0.0 "`γ` has to be greater than 0"
        return new{T}(γ)
    end
end

struct SquaredMax{T<:AbstractFloat} <: MaxOperator
    γ::T
    function SquaredMax(γ::T=1.0) where {T<:AbstractFloat}
        @assert γ > 0.0 "`γ` has to be greater than 0"
        return new{T}(γ)
    end
end

# Maximization functions
# ----------------------

#maximum(mo::MaxOperator, x::Vector{<:Number}) = max(mo, x...)
maximum(mo::Max, x::Vector{<:Number}) = maximum(x)

function maximum(mo::LeakyMax{T}, x::Vector{<:Number}) where {T}
    p = mo.p
    return (one(T) - p) * maximum(x) + p * fin_mean(x)
end

function maximum(mo::EntropyMax, x::Vector{<:Number})
    γ = mo.γ
    return x ./ γ .|> exp |> sum |> log |> m -> γ * m
end

function maximum(mo::SquaredMax, x::Vector{<:Number})
    γ = mo.γ
    q = project_in_simplex(x ./ γ , one(γ))
    return fin_dot(q, x) - (γ/2) * norm(q)^2
end

# Chain Rules
# -----------

function frule(::typeof(maximum), ::Max, x::Vector{T}) where {T<:Number}
    i = argmax(x)
    m = x[i]
    q = zeros(T, length(x))
    q[i] = one(T)
    return m, q
end

function frule(::typeof(maximum), mo::LeakyMax, x::Vector{<:Number})
    p = mo.p
    i = argmax(x)
    pcompl = (one(p) - p)
    m = pcompl * x[i] + p * fin_mean(x)
    n = count(el-> el > -Inf, x)
    q = Vector{typeof(p)}(undef, length(x))
    @. q = p / n * (x > -Inf)
    q[i] += pcompl
    return m, q
end

function frule(::typeof(maximum), mo::EntropyMax, x::Vector{<:Number})
    γ = mo.γ
    m = maximum(mo, x)
    q = exp.((x .-fin_mean(x)) ./ γ)
    q ./= sum(q)
    return m, q
end

function frule(::typeof(maximum), mo::SquaredMax, x::Vector{<:Number})
    γ = mo.γ
    q = project_in_simplex(x ./ γ , one(γ))
    m = fin_dot(q, x) - (γ/2) * norm(q)^2
    return m, q
end

"""
    max_argmax(mo::MaxOperator, x::Vector{<:Number})

Returns the maximum and the argmaxim (i.e., the gradient of the max) of 
a vector `x` using a given `MaxOperator`. 
"""
max_argmax(mo::MaxOperator, x::Vector{<:Number}) = frule(maximum, mo, x)


# Minimum
# -------

minimum(mo::Max, x::Vector{<:Number}) = minimum(x)
minimum(mo::MaxOperator, x::Vector{<:Number}) = -maximum(mo, -x)

function frule(::typeof(minimum), mo::MaxOperator, x::Vector{<:Number})
    m, q = frule(maximum, mo, -x)
    return -m, q
end

min_argmin(mo::MaxOperator, x::Vector{<:Number}) = frule(minimum, mo, x)

#=
min(mo::Max, x::Number...) = min(x...)

function min(mo::LeakyMax, x::Number...)
    !all(isfinite.(x)) && return min(x...)
    p = mo.p
    n = length(x)
    return (1.0 - p) * min(x...) + p * mean(x)
end

min(mo::EntropyMax, x::Number...) = ((x) ./ -mo.γ) .|> exp |> sum |> x -> -mo.γ * log(x)

function max_argmax!(::Max, x::VecOrMat{T}) where {T}
    i = argmax(x)
    m = x[i]
    x .= zero(T)
    x[i] = one(T)
    return m, x
end

function max_argmax!(mo::LeakyMax, x::VecOrMat{T}) where {T}
    # if there are infinities, just use regular max
    !all(isfinite.(x)) && return max_argmax!(Max(), x)
    p = mo.p
    n = length(x)
    i = argmax(x)
    m = (1.0 - p) * x[i] + p * mean(x)
    x .= p / n
    x[i] += (1.0 - p)
    return m, x
end

function max_argmax!(mo::EntropyMax, x::VecOrMat{T}) where {T}
    # if there are infinities, just use regular max
    !all(isfinite.(x)) && return max_argmax!(Max(), x)
    γ = mo.γ
    xm = mean(x)
    # substract mean, scale and exponent
    @. x = exp((x - xm) / γ)
    se = sum(x)
    x ./= se
    # FIXME: you cannot substract mean for log sum max
    return γ * log(se),  x
end

function max_argmax!(mo::SquaredMax, x::Vector)
    # if there are infinities, just use regular max
    !all(isfinite.(x)) && return max_argmax!(Max(), x)
    γ = mo.γ
    q = project_in_simplex(x ./ γ , 1.0)
    sqm = x ⋅ q - 0.5γ * norm(q)^2.0
    x .= q
    return sqm,  q
end

function max_argmax(mo::SquaredMax, x::Vector)
    γ = mo.γ
    q = project_in_simplex(x ./ γ , 1.0)
    sqm = x ⋅ q - 0.5γ * norm(q)^2.0
    return sqm,  q
end

max_argmax(mo::MaxOperator, x) = max_argmax!(mo, copy(x))

function min_argmin!(mo::MaxOperator, x)
    x .= -x
    m, am = max_argmax!(mo, x)
    return -m, am
end

min_argmin(mo::MaxOperator, x) = min_argmin!(mo, copy(x))
=#