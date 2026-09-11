using PoissonRandom, BenchmarkTools
using StableRNGs

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

# =============================================================================
# pois_rand across λ regimes
# =============================================================================

SUITE["pois_rand"] = BenchmarkGroup()

SUITE["pois_rand"]["small_lambda"] = @benchmarkable sum(
    pois_rand($rng, 2.0) for _ in 1:1000
)
SUITE["pois_rand"]["medium_lambda"] = @benchmarkable sum(
    pois_rand($rng, 30.0) for _ in 1:1000
)
SUITE["pois_rand"]["large_lambda"] = @benchmarkable sum(
    pois_rand($rng, 1.0e4) for _ in 1:1000
)
SUITE["pois_rand"]["huge_lambda"] = @benchmarkable sum(
    pois_rand($rng, 1.0e10) for _ in 1:1000
)

# PassthroughRNG (deterministic counting RNG, hot loop)
SUITE["pois_rand"]["passthrough"] = @benchmarkable sum(
    pois_rand(PassthroughRNG(), 5.0) for _ in 1:1000
)
