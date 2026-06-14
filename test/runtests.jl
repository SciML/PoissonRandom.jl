using PoissonRandom
import Distributions
using Test, Statistics
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "QA"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
    Pkg.instantiate()
    include(joinpath(@__DIR__, "qa", "qa.jl"))
else
    @safetestset "Sampler distribution tests" begin
        include("sampler_tests.jl")
    end

    @safetestset "BigFloat support" begin
        include("bigfloat_tests.jl")
    end

    @safetestset "PassthroughRNG dispatch" begin
        include("passthrough_rng_tests.jl")
    end

    @safetestset "Allocation Tests" begin
        include("alloc_tests.jl")
    end
end
