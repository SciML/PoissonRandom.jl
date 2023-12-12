using PoissonRandom, Aqua
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
