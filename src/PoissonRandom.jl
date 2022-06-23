module PoissonRandom

using Random

export pois_rand

count_rand(λ) = count_rand(Random.GLOBAL_RNG, λ)
function count_rand(rng::AbstractRNG, λ)
    n = 0
    c = randexp(rng)
    while c < λ
        n += 1
        c += randexp(rng)
    end
    return n
end

# Algorithm from:
#
#   J.H. Ahrens, U. Dieter (1982)
#   "Computer Generation of Poisson Deviates from Modified Normal Distributions"
#   ACM Transactions on Mathematical Software, 8(2):163-179
#
#   For μ sufficiently large, (i.e. >= 10.0)
#
ad_rand(λ) = ad_rand(Random.GLOBAL_RNG, λ)
function ad_rand(rng::AbstractRNG, λ)
    s = sqrt(λ)
    d = 6.0 * λ^2
    L = floor(Int, λ - 1.1484)
    # Step N
    G = λ + s * randn(rng)

    if G >= 0.0
        K = floor(Int, G)
        # Step I
        if K >= L
            return K
        end

        # Step S
        U = rand(rng)
        if d * U >= (λ - K)^3
            return K
        end

        # Step P
        px, py, fx, fy = procf(λ, K, s)

        # Step Q
        if fy * (1 - U) <= py * exp(px - fx)
            return K
        end
    end

    while true
        # Step E
        E = randexp(rng)
        U = 2.0 * rand(rng) - 1.0
        T = 1.8 + copysign(E, U)
        if T <= -0.6744
            continue
        end

        K = floor(Int, λ + s * T)
        px, py, fx, fy = procf(λ, K, s)
        c = 0.1069 / λ

        # Step H
        @fastmath if c * abs(U) <= py * exp(px + E) - fy * exp(fx + E)
            return K
        end
    end
end

# log(1+x)-x
# accurate ~2ulps for -0.227 < x < 0.315
function log1pmx_kernel(x::Float64)
    r = x / (x + 2.0)
    t = r * r
    w = @evalpoly(t,
                  6.66666666666666667e-1, # 2/3
                  4.00000000000000000e-1, # 2/5
                  2.85714285714285714e-1, # 2/7
                  2.22222222222222222e-1, # 2/9
                  1.81818181818181818e-1, # 2/11
                  1.53846153846153846e-1, # 2/13
                  1.33333333333333333e-1, # 2/15
                  1.17647058823529412e-1) # 2/17
    hxsq = 0.5 * x * x
    r * (hxsq + w * t) - hxsq
end

# use naive calculation or range reduction outside kernel range.
# accurate ~2ulps for all x
function log1pmx(x::Float64)
    if !(-0.7 < x < 0.9)
        return log1p(x) - x
    elseif x > 0.315
        u = (x - 0.5) / 1.5
        return log1pmx_kernel(u) - 9.45348918918356180e-2 - 0.5 * u
    elseif x > -0.227
        return log1pmx_kernel(x)
    elseif x > -0.4
        u = (x + 0.25) / 0.75
        return log1pmx_kernel(u) - 3.76820724517809274e-2 + 0.25 * u
    elseif x > -0.6
        u = (x + 0.5) * 2.0
        return log1pmx_kernel(u) - 1.93147180559945309e-1 + 0.5 * u
    else
        u = (x + 0.625) / 0.375
        return log1pmx_kernel(u) - 3.55829253011726237e-1 + 0.625 * u
    end
end

# Procedure F
function procf(λ, K::Int, s::Float64)
    # can be pre-computed, but does not seem to affect performance
    ω = 0.3989422804014327 / s
    b1 = 0.041666666666666664 / λ
    b2 = 0.3 * b1 * b1
    c3 = 0.14285714285714285 * b1 * b2
    c2 = b2 - 15.0 * c3
    c1 = b1 - 6.0 * b2 + 45.0 * c3
    c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3

    if K < 10
        px = -float(λ)
        py = λ^K / factorial(K)
    else
        δ = 0.08333333333333333 / K
        δ -= 4.8 * δ^3
        V = (λ - K) / K
        px = K * log1pmx(V) - δ # avoids need for table
        py = 0.3989422804014327 / sqrt(K)
    end
    X = (K - λ + 0.5) / s
    X2 = X^2
    fx = -0.5 * X2 # missing negation in pseudo-algorithm, but appears in fortran code.
    fy = ω * (((c3 * X2 + c2) * X2 + c1) * X2 + c0)
    return px, py, fx, fy
end

"""
```julia
pois_rand(λ)
pois_rand(rng::AbstractRNG, λ)
```

Generates Poisson(λ) distributed random numbers using a fast polyalgorithm.

## Examples

```julia
# Simple Poisson random
pois_rand(λ)

# Using another RNG
using RandomNumbers
rng = Xorshifts.Xoroshiro128Plus()
pois_rand(rng,λ)
```
"""
pois_rand(λ) = pois_rand(Random.GLOBAL_RNG, λ)
pois_rand(rng::AbstractRNG, λ) = λ < 6 ? count_rand(rng, λ) : ad_rand(rng, λ)

end # module
