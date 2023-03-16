# PoissonRandom.jl: Fast Poisson Random Numbers

PoissonRandom.jl is a component of the SciML ecosystem which allows
for fast generation of Poisson random numbers.

## Installation

To install PoissonRandom.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("PoissonRandom")
```

## Usage

```julia
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

## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to SciML.
- There are a few community forums:
    - the #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
    - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
    - on the [Julia Discourse forums](https://discourse.julialang.org)
    - see also [SciML Community page](https://sciml.ai/community/)

## Reproducibility
```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```
```@example
using Pkg # hide
Pkg.status() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>and using this machine and Julia version.</summary>
```
```@example
using InteractiveUtils # hide
versioninfo() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```
```@example
using Pkg # hide
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```
```@raw html
</details>
```
```@raw html
You can also download the 
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```
```@raw html
">project</a> file.
```
