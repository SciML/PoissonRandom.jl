using AllocCheck
using PoissonRandom
using Random

@testset "AllocCheck - Zero Allocations" begin
    @testset "count_rand" begin
        @check_allocs function test_count_rand(rng::TaskLocalRNG, λ::Float64)
            PoissonRandom.count_rand(rng, λ)
        end
        rng = Random.default_rng()
        @test test_count_rand(rng, 2.0) isa Int
        @test test_count_rand(rng, 5.0) isa Int
    end

    @testset "ad_rand" begin
        @check_allocs function test_ad_rand(rng::TaskLocalRNG, λ::Float64)
            PoissonRandom.ad_rand(rng, λ)
        end
        rng = Random.default_rng()
        @test test_ad_rand(rng, 10.0) isa Int
        @test test_ad_rand(rng, 50.0) isa Int
    end

    @testset "pois_rand" begin
        @check_allocs function test_pois_rand(rng::TaskLocalRNG, λ::Float64)
            pois_rand(rng, λ)
        end
        rng = Random.default_rng()
        @test test_pois_rand(rng, 2.0) isa Int
        @test test_pois_rand(rng, 10.0) isa Int
    end

    @testset "pois_rand with PassthroughRNG" begin
        @check_allocs function test_pois_rand_passthrough(rng::PassthroughRNG, λ::Float64)
            pois_rand(rng, λ)
        end
        passthrough = PassthroughRNG()
        @test test_pois_rand_passthrough(passthrough, 2.0) isa Int
        @test test_pois_rand_passthrough(passthrough, 10.0) isa Int
    end

    @testset "procf" begin
        @check_allocs function test_procf(λ::Float64, K::Int, s::Float64)
            PoissonRandom.procf(λ, K, s)
        end
        @test test_procf(10.0, 5, 3.162) isa NTuple{4,Float64}
        @test test_procf(10.0, 15, 3.162) isa NTuple{4,Float64}
    end
end
