# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.shrink.R
#
# Shrinkage estimators of entropy, mutual information and related quantities.

export getentropyshrinkage

"""
"""
function getlambda(n, normalizedcounts, target)

	# Unbiased estimator of variance of u
	varu = normalizedcounts .* (1 - normalizedcounts) / (n - 1)
	msp = sum((normalizedcounts - target).^2) # misspecification ???

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

lambdaFreqs
"""
function getfrequenciesshrinkage(counts, lambdaFreqs)

	target = 1 / length(counts) # Target is uniform distribution
	n = sum(counts)
	normalizedcounts = counts / n

	if lambdaFreqs == false
		lambdaFreqs = n == 1 || n == 0 ? 1 : getlambda(n, normalizedcounts, target)
	end

	return lambdaFreqs * target + (1 - lambdaFreqs) * normalizedcounts

end

"""
Calculates the shrinkage estimate for the entropy from the
observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

lambdaFreqs -

base - Int - The base of the logarithm, i.e. the units.
"""
function getentropyshrinkage(counts::Array{Float64,2}, base=2, lambdaFreqs=false)
	frequencies = getfrequenciesshrinkage(getfrequencies(counts), lambdaFreqs)
	return applyentropyformula(frequencies, base)
end
