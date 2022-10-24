# PoissonRandom.jl: Fast Poisson Random Numbers

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/PoissonRandom/stable/)

[![codecov](https://codecov.io/gh/SciML/PoissonRandom.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/PoissonRandom.jl)
[![Build Status](https://github.com/SciML/PoissonRandom.jl/workflows/CI/badge.svg)](https://github.com/SciML/PoissonRandom.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://docs.sciml.ai/PoissonRandom/stable/). Use the
[in-development documentation](https://docs.sciml.ai/PoissonRandom/dev/) for the version of
the documentation, which contains the unreleased features.

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
