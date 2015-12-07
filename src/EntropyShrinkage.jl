# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.shrink.R
#
# Shrinkage estimators of entropy, mutual information and related quantities.

"""
Calculates the shrinkage estimate for the entropy from the
observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

lambdaFreqs -

base - Int - The base of the logarithm, i.e. the units.

verbose -
"""
function entropyshrinkage(counts::Array{Float64,2}, lambdaFreqs=false, base=2, verbose=true)
	frequencies = frequenciesshrinkage(counts, lambdaFreqs, verbose)
	h = entropyPlugin(f, base)
end

function freqsShrink(y, lambdaFreqs, verbose)

	target = 1 / length(y) # Target is uniform distribution
	n = sum(y)
	u = y / n

	if !lambdaFreqs
		if n == 1 || n == 0
			lambdaFreqs = 1
		else
			lambdaFreqs = getLambdaShrink(n, u, target, verbose)
		end
	else
		if verbose
			# println("Specified shrinkage intensity: ", round(lambda, 4), "\n")
		end
	end

	uShrink = lambdaFreqs * target + (1 - lambdaFreqs) * u

end