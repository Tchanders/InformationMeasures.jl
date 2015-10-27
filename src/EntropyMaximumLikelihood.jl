# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.empirical.R
#
# Empirical estimators of entropy, mutual information and related quantities.

export entropymaximumlikelihood

function frequenciesmaximumlikelihood(frequencies)
	return convert(Array{Float64}, frequencies / sum(frequencies))
end

function entropymaximumlikelihood(counts, base=2)
	frequencies = frequenciesmaximumlikelihood(discretizecounts(counts))
	return entropyequation(frequencies, base)
end
