# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.empirical.R
#
# Empirical estimators of entropy, mutual information and related quantities.

export getentropymaximumlikelihood

"""
Calculates the maximum likelihood estimate for the true
probability distribution from the bin frequencies of the
observed values.

Parameters:

frequencies - Array{Int,2} - The bin frequencies of the observed
values.
"""
function getprobabilitiesmaximumlikelihood(frequencies::Array{Int,2})
	return frequencies / sum(frequencies)
end

"""
Calculates the maximum likelihood estimate for the entropy of a
set of observed values. First the observed values are converted
to frequencies of discrete bins, then the frequencies are
converted to probablities, then the probabilities are run
through the entropy formula.

Parameters:

values - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions of each value and n is the number of
values.

base - Int - The base of the logarithm, i.e. the units.
"""
function getentropymaximumlikelihood(values::Array{Float64,2}, base=2)
	probabilities = getprobabilitiesmaximumlikelihood(getfrequencies(values))
	return applyentropyformula(probabilities, base)
end
