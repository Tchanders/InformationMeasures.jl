# TODO: Add tests for joint frequencies

using EntropyEstimators
using Base.Test

d = 1
n = 100
arr = rand(d, n)
b = sqrt(n) # getfrequencies function uses sqrt(n) bins

include("testEntropyKnn.jl")
include("testEntropyMaximumLikelihood.jl")
include("testEntropyMillerMadow.jl")
include("testEntropyDirichlet.jl")
include("testEntropyShrinkage.jl")
