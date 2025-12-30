using PoissonRandom, Aqua, JET, ExplicitImports
using Random

@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(PoissonRandom)
    Aqua.test_ambiguities(PoissonRandom, recursive = false)
    Aqua.test_deps_compat(PoissonRandom)
    Aqua.test_piracies(PoissonRandom,
        treat_as_own = [])
    Aqua.test_project_extras(PoissonRandom)
    Aqua.test_stale_deps(PoissonRandom)
    Aqua.test_unbound_args(PoissonRandom)
    Aqua.test_undefined_exports(PoissonRandom)
end

@testset "JET static analysis" begin
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
end

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(PoissonRandom) === nothing
    @test check_no_stale_explicit_imports(PoissonRandom) === nothing
end
