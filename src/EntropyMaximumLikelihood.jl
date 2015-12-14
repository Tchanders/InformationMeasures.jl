# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.empirical.R
#
# Empirical estimators of entropy, mutual information and related quantities.

export getentropymaximumlikelihood

"""
Calculates the maximum likelihood estimate for the true
probability distribution from the observed counts.

Parameters:

counts - Array{Float64,1} - The observed bin frequencies.
"""
function getprobabilitiesmaximumlikelihood(counts::Array{Float64})
	return counts / sum(counts)
end

"""
Calculates the maximum likelihood estimate for the entropy from
the observed counts.

Parameters:

counts - dxn Array{Float64,2} - The observed counts.

base - Int - The base of the logarithm, i.e. the units.
"""
function getentropymaximumlikelihood(counts::Array{Float64,2}, base=2)
	probabilities = getprobabilitiesmaximumlikelihood(getfrequencies(counts))
	return applyentropyformula(probabilities, base)
end

"""
"""
function getmutualinformationmaximumlikelihood(counts1::Array{Float64,2}, counts2::Array{Float64,2}, base=2)
	jointfrequencies = getprobabilitiesmaximumlikelihood(
		makejointfrequencies(getfrequencies(counts1), getfrequencies(counts2))
	)
	

end
