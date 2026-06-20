using AllocCheck
using PoissonRandom
using Random
using Test

# AllocCheck's static LLVM-IR analysis is explicitly documented as not guaranteed
# stable across Julia versions, and on Julia prereleases the in-flux codegen makes
# it emit architecture-dependent false positives: on 1.13.0-rc1 the aarch64 path
# flags even `procf` (pure scalar Float64 arithmetic returning an `NTuple{4,Float64}`,
# which cannot heap-allocate) as allocating, while x86_64 reports zero. The runtime
# `@allocated == 0` check below is the ground-truth regression guard and is exercised
# on every version/architecture; the static `@check_allocs` guard is additionally run
# only on released Julia, where its result is reliable.
const RELIABLE_STATIC_ALLOC_CHECK = isempty(VERSION.prerelease)

# Assert `f(args...)` performs zero heap allocations at runtime (the real guard), and
# additionally pass it through AllocCheck's static analysis on released Julia.
macro test_zero_allocs(call)
    @assert call.head == :call
    f = call.args[1]
    args = call.args[2:end]
    return quote
        local f = $(esc(f))
        local args = ($(map(esc, args)...),)
        f(args...)                       # warm up / compile
        @test (@allocated f(args...)) == 0
        if RELIABLE_STATIC_ALLOC_CHECK
            local checked = AllocCheck.check_allocs(f, map(typeof, args))
            @test isempty(checked)
        end
    end
end

@testset "AllocCheck - Zero Allocations" begin
    rng = Random.default_rng()
    passthrough = PassthroughRNG()

    @testset "count_rand" begin
        @test_zero_allocs PoissonRandom.count_rand(rng, 2.0)
        @test_zero_allocs PoissonRandom.count_rand(rng, 5.0)
    end

    @testset "ad_rand" begin
        @test_zero_allocs PoissonRandom.ad_rand(rng, 10.0)
        @test_zero_allocs PoissonRandom.ad_rand(rng, 50.0)
    end

    @testset "pois_rand" begin
        @test_zero_allocs pois_rand(rng, 2.0)
        @test_zero_allocs pois_rand(rng, 10.0)
    end

    @testset "pois_rand with PassthroughRNG" begin
        @test_zero_allocs pois_rand(passthrough, 2.0)
        @test_zero_allocs pois_rand(passthrough, 10.0)
    end

    @testset "procf" begin
        @test_zero_allocs PoissonRandom.procf(10.0, 5, 3.162)
        @test_zero_allocs PoissonRandom.procf(10.0, 15, 3.162)
    end
end
