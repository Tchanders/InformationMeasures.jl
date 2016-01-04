# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.shrink.R
#
# Shrinkage estimators of entropy, mutual information and related quantities.

export getentropyshrinkage

"""
Calculates the shrinkage intensity, lambda

Parameters:

n - Int - The number of observed values.

normalizedvalues - The observed values adjusted so they sum to
1. This is equivalent to the maximum likelihood probabilities.

target - WIP
"""
function getlambda(n::Int, normalizedvalues, target)

	# Unbiased estimator of variance of u
	varu = normalizedvalues .* (1 - normalizedvalues) / (n - 1)
	msp = sum((normalizedvalues - target).^2) # misspecification ???

	# Estimate shrinkage intensity
	lambda = msp == 0 ? 1 : sum(varu) / msp
	
	# Make lambda be between 0 and 1 inclusive
	return lambda > 1 ? 1 : (lambda < 0 ? 0 : lambda)

end

"""
Calculates the shrinkage estimate for the true probability
distribution from the observed values.

Parameters:

values - Array{Float64} - The observed values.

lambda - WIP
"""
function getprobabilitiesshrinkage(values::Array{Float64,1})

	target = 1 / length(values) # Target is uniform distribution
	n = sum(values)
	normalizedvalues = values / n
	lambda = n == 1 || n == 0 ? 1 : getlambda(n, normalizedvalues, target)

	return lambda * target + (1 - lambda) * normalizedvalues

end
function getprobabilitiesshrinkage(values::Array{Float64,1}, lambda::Number)

	target = 1 / length(values) # Target is uniform distribution
	n = sum(values)
	normalizedcounts = values / n

	return lambda * target + (1 - lambda) * normalizedcounts

end

"""
Calculates the shrinkage estimate for the entropy of a set of
observed values. First the observed values are converted to
frequencies of discrete bins, then the frequencies are converted
to probablities using shrinkage estimation, then the probabilities
are run through the entropy formula.

Parameters:

values - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions of each value and n is the number of
values.

base - Int - The base of the logarithm, i.e. the units.

lambda - WIP
"""
function getentropyshrinkage(counts::Array{Float64,2}, base=2, lambda=false)
	probabilities = getprobabilitiesshrinkage(getfrequencies(counts), lambda)
	return applyentropyformula(probabilities, base)
end
