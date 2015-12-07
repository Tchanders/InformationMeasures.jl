# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.empirical.R
#
# Empirical estimators of entropy, mutual information and related quantities.

export entropymaximumlikelihood

"""
Calculates the maximum likelihood estimate for the true
probability distribution from the observed counts.

Parameters:

counts - Array{Float64,1} - The observed bin frequencies.
"""
function frequenciesmaximumlikelihood(counts::Array{Float64,1})
	return convert(Array{Float64}, counts / sum(counts))
end

"""
Calculates the maximum likelihood estimate for the entropy from
the observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

base - Int - The base of the logarithm, i.e. the units.
"""
function entropymaximumlikelihood(counts::Array{Float64,2}, base=2)
	frequencies = frequenciesmaximumlikelihood(discretizecounts(counts))
	return entropyequation(frequencies, base)
end
