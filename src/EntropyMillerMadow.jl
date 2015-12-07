# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.MillerMadow.R
#
# Miller-Madow entropy estimator (1955)

export entropymillermadow

"""
Calculates the Miller Madow estimate for the entropy from the
observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

base - Int - The base of the logarithm, i.e. the units.
"""
function entropymillermadow(counts::Array{Float64,2}, base=2)
	frequencies = frequenciesmaximumlikelihood(discretizecounts(counts))
	constant = (countnz(frequencies) - 1) / (2 * length(frequencies))
	return entropyequation(frequencies, base) + constant
end
