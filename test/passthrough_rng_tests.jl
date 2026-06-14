using PoissonRandom
using Random: Random, UInt52Raw
using Test

prng = PassthroughRNG()
# The CUDA.jl @device_override Random.randexp(::AbstractRNG) shadows our
# specific Random.randexp(::PassthroughRNG) on the GPU because Julia's
# OverlayMethodTable returns overlay matches without consulting the base
# table when the overlay fully covers the signature. The override body
# then calls these against PassthroughRNG; if they MethodError, kernel
# compilation fails with InvalidIRError on jl_f_throw_methoderror.
@test Random.rng_native_52(prng) === UInt64
@test Random.rand(prng, UInt52Raw()) isa UInt64
@test Random.rand(prng, UInt64) isa UInt64
@test Random.rand(prng, Float32) isa Float32
@test Random.rand(prng, Float64) isa Float64
