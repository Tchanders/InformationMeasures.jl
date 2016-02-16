# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.MillerMadow.R
#
# Miller-Madow entropy estimator (1955)

export getentropymillermadow

"""
Calculates the Miller Madow estimate for the entropy of a set
of observed values. First the observed values are converted
to frequencies of discrete bins, then the frequencies are
converted to probablities (using maximum likelihood), then
the probabilities are run through the entropy formula, then
a constant is added to correct for bias.

Parameters:

values - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions of each value and n is the number of
values.

base - Int - The base of the logarithm, i.e. the units.
"""
function getentropymillermadow(values::Array{Float64,2}, base=2)
	probabilities = getprobabilitiesmaximumlikelihood(getfrequencies(values))
	constant = (countnz(probabilities) - 1) / (2 * length(probabilities))
	return applyentropyformula(probabilities, base) + constant
end
function getentropymillermadow(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, base=2)
	probabilities = getprobabilitiesmaximumlikelihood(getjointfrequencies(valuesX, valuesY))
	constant = (countnz(probabilities) - 1) / (2 * length(probabilities))
	return applyentropyformula(probabilities, base)
end
