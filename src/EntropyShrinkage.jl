# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.shrink.R
#
# Shrinkage estimators of entropy, mutual information and related quantities.

export getentropyshrinkage

"""
Calculates the shrinkage intensity, lambda

Parameters:

n - The total number of observations

counts -

target -
"""
function getlambda(n, normalizedcounts, target)

	# Unbiased estimator of variance of u
	varu = counts .* (1 - counts) / (n - 1)
	msp = sum((counts - target).^2) # misspecification ???

	# Estimate shrinkage intensity
	lambda = msp == 0 ? 1 : sum(varu) / msp
	
	# Make lambda be between 0 and 1 inclusive
	return lambda > 1 ? 1 : (lambda < 0 ? 0 : lambda)

end

"""
Calculates an estimate for the true probability distribution
from the observed counts.

Parameters:

counts - Array{Float64,1} - The observed bin frequencies.

lambda -
"""
function getprobabilitiesshrinkage(counts::Array{Float64,1})

	target = 1 / length(counts) # Target is uniform distribution
	n = sum(counts)
	normalizedcounts = counts / n
	lambda = n == 1 || n == 0 ? 1 : getlambda(n, normalizedcounts, target)

	return lambda * target + (1 - lambda) * normalizedcounts

end
function getprobabilitiesshrinkage(counts::Array{Float64,1}, lambda::Number)

	target = 1 / length(counts) # Target is uniform distribution
	n = sum(counts)
	normalizedcounts = counts / n

	return lambda * target + (1 - lambda) * normalizedcounts

end

"""
Calculates the shrinkage estimate for the entropy from the
observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

base - Int - The base of the logarithm, i.e. the units.

lambda -
"""
function getentropyshrinkage(counts::Array{Float64,2}, base=2, lambdaFreqs=false)
	probabilities = getprobabilitiesshrinkage(getfrequencies(counts), lambdaFreqs)
	return applyentropyformula(probabilities, base)
end
