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

export getentropydirichlet

"""
Calculates the Dirichlet prior Bayesian estimate for the true
probability distribution from the bin frequencies of the
observed values.

Parameters:

frequencies - Array{Int,1} - The observed values.

a - Float64 - The Dirichlet prior.
"""
function getprobabilitiesdirichlet(frequencies::Array{Int,1}, a::Number)
	a = fill(a, length(frequencies))
	return (frequencies + a) / (sum(frequencies) + sum(a))
end

"""
Calculates the Dirichlet prior Bayesian estimate for the entropy
of a set of observed values. First the observed values are
converted to frequencies of discrete bins, then the frequencies
are converted to probablities using a Dirichlet prior, then the
probabilities are run through the entropy formula.

Parameters:

values - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions of each value and n is the number of
values.

a - Any - The Dirichlet prior.

base - Int - The base of the logarithm, i.e. the units.
"""
function getentropydirichlet(values::Array{Float64,2}, a::Any, base=2)
	probabilities = getprobabilitiesdirichlet(getfrequencies(values), a)
	return applyentropyformula(probabilities, base)
end
