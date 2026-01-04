module PoissonRandom

using Random: Random, AbstractRNG, randexp
using LogExpFunctions: log1pmx
using PrecompileTools: @compile_workload

export pois_rand, PassthroughRNG

# GPU-compatible Poisson sampling via PassthroughRNG
struct PassthroughRNG <: AbstractRNG end

Random.rand(rng::PassthroughRNG) = rand()
Random.randexp(rng::PassthroughRNG) = randexp()
Random.randn(rng::PassthroughRNG) = randn()

count_rand(λ::Real) = count_rand(Random.GLOBAL_RNG, λ)
function count_rand(rng::AbstractRNG, λ::Real)
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
ad_rand(λ::Real) = ad_rand(Random.GLOBAL_RNG, λ)
function ad_rand(rng::AbstractRNG, λ::Real)
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
    return
end

# Procedure F
function procf(λ::Real, K::Int, s::Real)
    # can be pre-computed, but does not seem to affect performance
    INV_SQRT_2PI = inv(sqrt(2pi))
    ω = INV_SQRT_2PI / s
    b1 = inv(24) / λ
    b2 = 0.3 * b1 * b1
    c3 = inv(7) * b1 * b2
    c2 = b2 - 15 * c3
    c1 = b1 - 6 * b2 + 45 * c3
    c0 = 1 - b1 + 3 * b2 - 15 * c3

    if K < 10
        px = -float(λ)
        py = λ^K / prod(2:K)
    else
        δ = inv(12) / K
        δ -= 4.8 * δ^3
        V = (λ - K) / K
        px = K * log1pmx(V) - δ # avoids need for table
        py = INV_SQRT_2PI / sqrt(K)
    end
    X = (K - λ + 0.5) / s
    X2 = X^2
    fx = X2 / -2 # missing negation in pseudo-algorithm, but appears in fortran code.
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
pois_rand(rng, λ)

# Simple Poisson random on GPU
pois_rand(PoissonRandom.PassthroughRNG(), λ)
```
"""
pois_rand(λ::Real) = pois_rand(Random.GLOBAL_RNG, λ)
pois_rand(rng::AbstractRNG, λ::Real) = λ < 6 ? count_rand(rng, λ) : ad_rand(rng, λ)

@compile_workload begin
    # Precompile the most common code paths
    # Small λ uses count_rand, large λ uses ad_rand
    pois_rand(3.0)   # count_rand path (λ < 6)
    pois_rand(50.0)  # ad_rand path (λ >= 6)
    # PassthroughRNG for GPU compatibility
    pois_rand(PassthroughRNG(), 3.0)
    pois_rand(PassthroughRNG(), 50.0)
end

end # module
