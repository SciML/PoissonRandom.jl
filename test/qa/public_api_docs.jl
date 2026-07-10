using PoissonRandom
using Test

@testset "public API documentation" begin
    public_names = filter(!=(:PoissonRandom), names(PoissonRandom; all = false, imported = false))
    @test Set(public_names) == Set([:PassthroughRNG, :pois_rand])

    for name in public_names
        binding = Docs.Binding(PoissonRandom, name)
        @test Docs.hasdoc(binding)
    end

    api_page = read(joinpath(pkgdir(PoissonRandom), "docs", "src", "pois_rand.md"), String)
    for name in public_names
        @test occursin(string(name), api_page)
    end
end
