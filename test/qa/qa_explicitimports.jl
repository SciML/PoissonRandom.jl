using PoissonRandom, ExplicitImports
using Test

@test check_no_implicit_imports(PoissonRandom) === nothing
@test check_no_stale_explicit_imports(PoissonRandom) === nothing
