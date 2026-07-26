using SciMLTesting, PoissonRandom, JET

run_qa(
    PoissonRandom;
    ei_kwargs = (;
        # `rng_native_52` is the stdlib sampler hook required by CUDA's device
        # `randexp` overlay for `PassthroughRNG`; Random has no public replacement.
        all_qualified_accesses_are_public = (; ignore = (:rng_native_52,)),
    ),
)
