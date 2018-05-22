# PoissonRandom

[![Build Status](https://travis-ci.org/ChrisRackauckas/PoissonRandom.jl.svg?branch=master)](https://travis-ci.org/ChrisRackauckas/PoissonRandom.jl)

[![Coverage Status](https://coveralls.io/repos/ChrisRackauckas/PoissonRandom.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ChrisRackauckas/PoissonRandom.jl?branch=master)

[![codecov.io](http://codecov.io/github/ChrisRackauckas/PoissonRandom.jl/coverage.svg?branch=master)](http://codecov.io/github/ChrisRackauckas/PoissonRandom.jl?branch=master)

## Benchmark

```julia
using RandomNumbers, Distributions, BenchmarkTools, StaticArrays,
      RecursiveArrayTools, Plots, PoissonRandom
labels = ["count_rand","ad_rand","pois_rand","Distributions.jl"]
rng = Xorshifts.Xoroshiro128Plus()

function n_count(λ,rng,n)
  tmp = 0
  for i in 1:n
    tmp += PoissonRandom.count_rand(λ,rng)
  end
end

function n_pois(λ,rng,n)
  tmp = 0
  for i in 1:n
    tmp += pois_rand(λ,rng)
  end
end

function n_ad(λ,rng,n)
  tmp = 0
  for i in 1:n
    tmp += PoissonRandom.ad_rand(λ,rng)
  end
end

function n_dist(λ,rng,n)
  tmp = 0
  for i in 1:n
    tmp += rand(Poisson(λ))
  end
end

function time_λ(λ,rng,n)
  t1 = @elapsed n_count(λ,rng,n)
  t2 = @elapsed n_ad(λ,rng,n)
  t3 = @elapsed n_pois(λ,rng,n)
  t4 = @elapsed n_dist(λ,rng,n)
  @SArray [t1,t2,t3,t4]
end

# Compile
time_λ(5,rng,5000000)
# Run with a bunch of λ
times = VectorOfArray([time_λ(n,rng,5000000) for n in 1:20])'
plot(times,labels = labels, lw = 3)
```
