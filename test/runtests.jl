using PoissonRandom
import Distributions
using Test, Statistics

n_tsamples = 10^5

function test_samples(rand_func,
                      distr::Distributions.DiscreteUnivariateDistribution,
                      n::Int;                                # number of samples to generate
                      q::Float64 = 1.0e-8,                     # confidence interval, 1 - q as confidence
                      verbose::Bool = false)                   # show intermediate info (for debugging)

    # The basic idea
    # ------------------
    #   Generate n samples, and count the occurences of each value within a reasonable range.
    #   For each distinct value, it computes an confidence interval of the counts
    #   and checks whether the count is within this interval.
    #
    #   If the distribution has a bounded range, it also checks whether
    #   the samples are all within this range.
    #
    #   By setting a small q, we ensure that failure of the tests rarely
    #   happen in practice.
    #
    λ = distr.λ
    n > 1 || error("The number of samples must be greater than 1.")
    0.0 < q < 0.1 || error("The value of q must be within the open interval (0.0, 0.1).")

    # determine the range of values to examine
    vmin = minimum(distr)
    vmax = maximum(distr)

    rmin = floor(Int, quantile(distr, 0.00001))::Int
    rmax = floor(Int, quantile(distr, 0.99999))::Int
    m = rmax - rmin + 1  # length of the range
    p0 = Distributions.pdf.((distr,), rmin:rmax)  # reference probability masses
    @assert length(p0) == m

    # determine confidence intervals for counts:
    # with probability q, the count will be out of this interval.
    #
    clb = Vector{Int}(undef, m)
    cub = Vector{Int}(undef, m)
    for i in 1:m
        bp = Distributions.Binomial(n, p0[i])
        clb[i] = floor(Int, quantile(bp, q / 2))
        cub[i] = ceil(Int, Distributions.cquantile(bp, q / 2))
        @assert cub[i] >= clb[i]
    end

    # generate samples
    samples = [rand_func(λ) for i in 1:n]
    @assert length(samples) == n

    # scan samples and get counts
    cnts = zeros(Int, m)
    for i in 1:n
        @inbounds si = samples[i]
        if rmin <= si <= rmax
            cnts[si - rmin + 1] += 1
        else
            vmin <= si <= vmax ||
                error("Sample value out of valid range.")
        end
    end

    # check the counts
    for i in 1:m
        verbose && println("v = $(rmin+i-1) ==> ($(clb[i]), $(cub[i])): $(cnts[i])")
        clb[i] <= cnts[i] <= cub[i] ||
            error("The counts are out of the confidence interval.")
    end
    return samples
end

println("testing count random sampler")
for λ in [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]
    test_samples(PoissonRandom.count_rand, Distributions.Poisson(λ), n_tsamples)
end

println("testing ad random sampler")
for λ in [5.0, 10.0, 15.0, 20.0, 30.0]
    test_samples(PoissonRandom.ad_rand, Distributions.Poisson(λ), n_tsamples)
end

println("testing mixed random sampler")
for λ in [5.0, 10.0, 15.0, 20.0, 30.0]
    test_samples(pois_rand, Distributions.Poisson(λ), n_tsamples)
end
