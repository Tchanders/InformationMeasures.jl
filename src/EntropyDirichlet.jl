# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.Dirichlet.R
#
# Dirichlet prior Bayesian estimator of entropy

# estimate entropy based on Dirichlet-multinomial pseudocount model

# y:  a vector of counts (may include zeros)
# a:  pseudocount per bin

# some choices for a:
# a = 0          :   empirical estimate
# a = 1          :   Laplace
# a = 1/2        :   Jeffreys
# a = 1/m        :   Schurmann-Grassberger  (m: number of bins)
# a = sqrt(n)/m  :   minimax

export entropydirichlet

"""
Calculates an estimate for the true probability distribution
from the observed frequencies.

Parameters:

frequencies - Array{Float64,1} - The observed bin frequencies.

a - Float64 - The Dirichlet prior.
"""
function frequenciesdirichlet(frequencies::Array{Float64,1}, a::Any)
	a = fill(a, length(frequencies))
	return (frequencies + a) / (sum(frequencies) + sum(a))
end

"""
Calculates the Dirichlet prior Bayesian estimate for the entropy
from the observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

a - Any - The Dirichlet prior.

base - Int - The base of the logarithm, i.e. the units.
"""
function entropydirichlet(counts::Array{Float64,2}, a::Any, base=2)
	frequencies = frequenciesdirichlet(discretizecounts(counts), a)
	return entropyequation(frequencies, base)
end
