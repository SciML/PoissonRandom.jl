module PoissonRandom

using Random
using LogExpFunctions: log1pmx
using SpecialFunctions: loggamma

export pois_rand, PassthroughRNG

# GPU-compatible Poisson sampling PassthroughRNG
struct PassthroughRNG <: AbstractRNG end

rand(rng::PassthroughRNG) = Random.rand()
randexp(rng::PassthroughRNG) = Random.randexp()
randn(rng::PassthroughRNG) = Random.randn()

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
    d = 6 * λ^2
    L = floor(Int, λ - 1.1484)
    # Step N
    G = λ + s * randn(rng)

    if G >= 0
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
        U = 2 * rand(rng) - 1
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

# Procedure F
function procf(λ, K::Int, s::Float64)
    INV_SQRT_2PI = 0.3989422804014327  # 1/sqrt(2π)
    ω = INV_SQRT_2PI / s
    b1 = 1 / (24 * λ)
    b2 = 0.3 * b1^2
    c3 = b1 * b2 / 7
    c2 = b2 - 15 * c3
    c1 = b1 - 6 * b2 + 45 * c3
    c0 = 1 - b1 + 3 * b2 - 15 * c3

    if K < 10
        px = -λ
        log_py = K * log(λ) - loggamma(K + 1)  # log(K!) via loggamma
        py = exp(log_py)
    else
        δ = 1 / (12 * K)
        δ -= 4.8 * δ^3
        V = (λ - K) / K
        px = K * log1pmx(V) - δ  # avoids need for table
        py = INV_SQRT_2PI / sqrt(K)
    end
    X = (K - λ + 0.5) / s
    X2 = X^2
    fx = -X2 / 2 # missing negation in pseudo-algorithm, but appears in fortran code.
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
# Simple Poisson random which works on GPU
pois_rand(λ)

# Using RNG
using RandomNumbers
rng = Xorshifts.Xoroshiro128Plus()
pois_rand(rng, λ)
```
"""
pois_rand(λ) = pois_rand(Random.GLOBAL_RNG, λ)
pois_rand(rng::AbstractRNG, λ) = λ < 6 ? count_rand(rng, λ) : ad_rand(rng, λ)

end # module
