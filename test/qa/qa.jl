using SciMLTesting, PoissonRandom, JET, Test

run_qa(
    PoissonRandom;
    explicit_imports = true,
    ei_kwargs = (;
        # default_rng / rng_native_52 are Random stdlib internals (not public in
        # Random): default_rng is the standard default-RNG accessor; rng_native_52 is
        # extended for the PassthroughRNG sampler chain.
        all_qualified_accesses_are_public = (; ignore = (:default_rng, :rng_native_52)),
    ),
)
