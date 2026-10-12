using PoissonRandom
using Random
using Statistics
using Test

# procf used to return Tuple{Union{Float32, Float64}, ...} for Float32 rates
# because its Float64 constants promoted only part of the tuple. Destructuring
# that Union-typed tuple in ad_rand forces a dynamic call, which cannot compile
# in GPU kernels (https://github.com/SciML/PoissonRandom.jl/issues/102).
@testset "procf type stability" begin
    for T in (Float32, Float64)
        λ = T(15)
        s = sqrt(λ)
        @test only(Base.return_types(PoissonRandom.procf, (T, Int, T))) === NTuple{4, T}
        @test @inferred(PoissonRandom.procf(λ, 5, s)) isa NTuple{4, T}
        @test @inferred(PoissonRandom.procf(λ, 15, s)) isa NTuple{4, T}
    end

    # No Union-typed value may appear in ad_rand for Float32 rates: each one is
    # a dynamic call site in a GPU kernel.
    ir = first(
        Base.code_typed(
            PoissonRandom.ad_rand, (Xoshiro, Float32);
            optimize = false
        )
    )[1]
    @test !any(t -> t isa Union, ir.slottypes)
    @test !any(t -> t isa Union, ir.ssavaluetypes)
end

@testset "Float32 sampling" begin
    rng = Xoshiro(1234)
    λ = 15.0f0
    samples = [pois_rand(rng, λ) for _ in 1:10_000]
    @test all(s -> s isa Integer && s >= 0, samples)
    # Poisson(λ) has mean = var = λ; standard errors here are ≈0.04 and ≈0.22.
    @test abs(mean(samples) - λ) < 0.5
    @test abs(var(samples) - λ) < 1.5

    @test pois_rand(PassthroughRNG(), λ) isa Integer
end
