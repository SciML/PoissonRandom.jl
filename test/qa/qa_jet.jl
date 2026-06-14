using PoissonRandom, JET
using Random
using Test

@testset "Type stability" begin
    JET.@test_opt target_modules = (PoissonRandom,) pois_rand(10.0)
    JET.@test_opt target_modules = (PoissonRandom,) pois_rand(Random.default_rng(), 10.0)
    JET.@test_opt target_modules = (PoissonRandom,) pois_rand(PassthroughRNG(), 10.0)
end

@testset "Error analysis" begin
    JET.@test_call target_modules = (PoissonRandom,) pois_rand(10.0)
    JET.@test_call target_modules = (PoissonRandom,) pois_rand(Random.default_rng(), 10.0)
    JET.@test_call target_modules = (PoissonRandom,) pois_rand(PassthroughRNG(), 10.0)
end
