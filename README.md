# PoissonRandom

[![Build Status](https://travis-ci.org/JuliaDiffEq/PoissonRandom.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/PoissonRandom.jl)

## Usage

```julia
Pkg.add("PoissonRandom")

# Simple Poisson random
pois_rand(λ)

# Using another RNG
using RandomNumbers
rng = Xorshifts.Xoroshiro128Plus()
pois_rand(rng,λ)
```

## Implementation

It mixes two methods. A simple count method and a method from a normal approximation.
See [this blog post for details](https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/).

## Benchmark

```julia
using RandomNumbers, Distributions, BenchmarkTools, StaticArrays,
      RecursiveArrayTools, Plots, PoissonRandom
labels = ["count_rand","ad_rand","pois_rand","Distributions.jl"]
rng = Xorshifts.Xoroshiro128Plus()

function n_count(rng,λ,n)
  tmp = 0
  for i in 1:n
    tmp += PoissonRandom.count_rand(rng,λ)
  end
end

function n_pois(rng,λ,n)
  tmp = 0
  for i in 1:n
    tmp += pois_rand(rng,λ)
  end
end

function n_ad(rng,λ,n)
  tmp = 0
  for i in 1:n
    tmp += PoissonRandom.ad_rand(rng,λ)
  end
end

function n_dist(λ,n)
  tmp = 0
  for i in 1:n
    tmp += rand(Poisson(λ))
  end
end

function time_λ(rng,λ,n)
  t1 = @elapsed n_count(rng,λ,n)
  t2 = @elapsed n_ad(rng,λ,n)
  t3 = @elapsed n_pois(rng,λ,n)
  t4 = @elapsed n_dist(λ,n)
  @SArray [t1,t2,t3,t4]
end

# Compile
time_λ(rng,5,5000000)
# Run with a bunch of λ
times = VectorOfArray([time_λ(rng,n,5000000) for n in 1:20])'
plot(times,labels = labels, lw = 3)
```

![benchmark result](https://user-images.githubusercontent.com/1814174/40387004-1e377776-5dc0-11e8-88a2-2d9cb12db984.png)

So this package ends up about 30% or so faster than Distributions.jl (the method
at the far edge is λ-independent so that goes on forever).
