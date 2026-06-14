using PoissonRandom
using Statistics
using Test

@testset "count_rand with BigFloat (λ < 6)" begin
    for _ in 1:100
        result = pois_rand(BigFloat(3.0))
        @test result isa Integer
        @test result >= 0
    end
end
@testset "ad_rand with BigFloat (λ >= 6)" begin
    for _ in 1:100
        result = pois_rand(BigFloat(15.0))
        @test result isa Integer
        @test result >= 0
    end
end
@testset "statistical validity with BigFloat" begin
    n = 10000
    λ = BigFloat(10.0)
    samples = [pois_rand(λ) for _ in 1:n]
    sample_mean = mean(samples)
    @test abs(sample_mean - Float64(λ)) < 3 * sqrt(Float64(λ))
end
