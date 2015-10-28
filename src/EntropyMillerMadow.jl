# Julia port of https://cran.r-project.org/web/packages/entropy/
#
# entropy.MillerMadow.R
#
# Miller-Madow entropy estimator (1955)

export entropymillermadow

function entropymillermadow(counts, base=2)
	frequencies = frequenciesmaximumlikelihood(discretizecounts(counts))
	c = (countnz(frequencies) - 1) / (2 * length(frequencies))
	return entropyequation(counts, base) + c
end
