var documenterSearchIndex = {"docs":
[{"location":"pois_rand/#Poisson-Random-API","page":"Poisson Random API","title":"Poisson Random API","text":"","category":"section"},{"location":"pois_rand/","page":"Poisson Random API","title":"Poisson Random API","text":"pois_rand","category":"page"},{"location":"pois_rand/#PoissonRandom.pois_rand","page":"Poisson Random API","title":"PoissonRandom.pois_rand","text":"pois_rand(λ)\npois_rand(rng::AbstractRNG, λ)\n\nGenerates Poisson(λ) distributed random numbers using a fast polyalgorithm.\n\nExamples\n\n# Simple Poisson random\npois_rand(λ)\n\n# Using another RNG\nusing RandomNumbers\nrng = Xorshifts.Xoroshiro128Plus()\npois_rand(rng, λ)\n\n\n\n\n\n","category":"function"},{"location":"#PoissonRandom.jl:-Fast-Poisson-Random-Numbers","page":"Home","title":"PoissonRandom.jl: Fast Poisson Random Numbers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PoissonRandom.jl is a component of the SciML ecosystem which allows for fast generation of Poisson random numbers.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install PoissonRandom.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"PoissonRandom\")","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"# Simple Poisson random\npois_rand(λ)\n\n# Using another RNG\nusing RandomNumbers\nrng = Xorshifts.Xoroshiro128Plus()\npois_rand(rng, λ)","category":"page"},{"location":"#Implementation","page":"Home","title":"Implementation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"It mixes two methods. A simple count method and a method from a normal approximation. See this blog post for details.","category":"page"},{"location":"#Benchmark","page":"Home","title":"Benchmark","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using RandomNumbers,\n      Distributions, BenchmarkTools, StaticArrays,\n      RecursiveArrayTools, Plots, PoissonRandom\nlabels = [\"count_rand\", \"ad_rand\", \"pois_rand\", \"Distributions.jl\"]\nrng = Xorshifts.Xoroshiro128Plus()\n\nfunction n_count(rng, λ, n)\n    tmp = 0\n    for i in 1:n\n        tmp += PoissonRandom.count_rand(rng, λ)\n    end\nend\n\nfunction n_pois(rng, λ, n)\n    tmp = 0\n    for i in 1:n\n        tmp += pois_rand(rng, λ)\n    end\nend\n\nfunction n_ad(rng, λ, n)\n    tmp = 0\n    for i in 1:n\n        tmp += PoissonRandom.ad_rand(rng, λ)\n    end\nend\n\nfunction n_dist(λ, n)\n    tmp = 0\n    for i in 1:n\n        tmp += rand(Poisson(λ))\n    end\nend\n\nfunction time_λ(rng, λ, n)\n    t1 = @elapsed n_count(rng, λ, n)\n    t2 = @elapsed n_ad(rng, λ, n)\n    t3 = @elapsed n_pois(rng, λ, n)\n    t4 = @elapsed n_dist(λ, n)\n    @SArray [t1, t2, t3, t4]\nend\n\n# Compile\ntime_λ(rng, 5, 5000000)\n# Run with a bunch of λ\ntimes = VectorOfArray([time_λ(rng, n, 5000000) for n in 1:20])'\nplot(times, labels = labels, lw = 3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: benchmark result)","category":"page"},{"location":"","page":"Home","title":"Home","text":"So this package ends up about 30% or so faster than Distributions.jl (the method at the far edge is λ-independent, so that goes on forever).","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please refer to the SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages for guidance on PRs, issues, and other matters relating to contributing to SciML.\nSee the SciML Style Guide for common coding practices and other style decisions.\nThere are a few community forums:\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Slack\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Zulip\nOn the Julia Discourse forums\nSee also SciML Community page","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status(; mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nusing Markdown\nversion = TOML.parse(read(\"../../Project.toml\", String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\", String))[\"name\"]\nlink_manifest = \"https://github.com/SciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n                \"/assets/Manifest.toml\"\nlink_project = \"https://github.com/SciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n               \"/assets/Project.toml\"\nMarkdown.parse(\"\"\"You can also download the\n[manifest]($link_manifest)\nfile and the\n[project]($link_project)\nfile.\n\"\"\")","category":"page"}]
}
