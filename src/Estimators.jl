# Functions for estimating probabilities using various estimators

# Estimators are described in:
# Hausser, Jean; Strimmer, Korbinian (2009-01-01).
# "Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks"
# https://arxiv.org/abs/0811.3579

# R implementation of estimators:
# https://cran.r-project.org/web/packages/entropy/

# Parameters:
# 	- frequencies, integer array
# 	- prior, number
function get_probabilities_dirichlet(frequencies, prior)
	prior_frequencies = fill(prior, size(frequencies))
	return (frequencies + prior_frequencies) / (sum(frequencies) + sum(prior_frequencies))
end

# Parameters:
# 	- frequencies, integer array
function get_probabilities_maximum_likelihood(frequencies)
	return frequencies / sum(frequencies)
end

# TODO: allow pre-calculated target to be passed in instead of uniform distribution
# Parameters:
# 	- frequencies, integer array
# 	- lambda, void
function get_probabilities_shrinkage(frequencies, lambda::Nothing)
	target = get_uniform_distribution(frequencies)
	n = sum(frequencies)
	normalized_frequencies = frequencies / n
	calculated_lambda = get_lambda(normalized_frequencies, target, n)
	return apply_shrinkage_formula(normalized_frequencies, target, calculated_lambda)
end
# Parameters:
# 	- frequencies, integer array
# 	- lambda, number
function get_probabilities_shrinkage(frequencies, lambda)
	target = get_uniform_distribution(frequencies)
	normalized_frequencies = get_normalized_frequencies(frequencies)
	return apply_shrinkage_formula(normalized_frequencies, target, lambda)
end

function apply_shrinkage_formula(normalized_frequencies, target, lambda)
	return lambda * target .+ (1 - lambda) * normalized_frequencies
end

function get_uniform_distribution(frequencies)
	return 1 / length(frequencies)
end

function get_normalized_frequencies(frequencies)
	return frequencies / sum(frequencies)
end

function get_lambda(normalized_frequencies, target, n)
	if n == 1 || n == 0
		return 1.0
	end
	# Add better comments about varu and msp
	# Unbiased estimator of variance of u
	varu = normalized_frequencies .* (1 .- normalized_frequencies) / (n - 1)
	msp = sum((normalized_frequencies .- target).^2) # misspecification ???

	# Estimate shrinkage intensity
	lambda = msp == 0 ? 1.0 : sum(varu) / msp

	# Make lambda be between 0 and 1 inclusive
	return lambda > 1 ? 1.0 : (lambda < 0 ? 0.0 : lambda)

end
function get_lambda(frequencies, get_target = get_uniform_distribution)
	n = sum(frequencies)
	if n == 1 || n == 0
		return 1
	end

	target = get_target(frequencies)
	normalized_frequencies = frequencies / n

	return get_lambda(normalized_frequencies, target, n)
end
