# The probability functions for the different estimators

# TODO: Reorganize shrinkage

# Parameters:
# 	Normal:
# 	- frequencies, integer array
# 	- a, number (currently optional but shouldn't be)
# 	Optional:
# 	Keyword:
function get_probabilities_dirichlet(frequencies, prior)
	prior = fill(prior, size(frequencies))
	return (frequencies + prior) / (sum(frequencies) + sum(prior))
end

# Parameters:
# 	Normal:
# 	- frequencies, integer array
# 	Optional:
# 	Keyword:
function get_probabilities_maximum_likelihood(frequencies)
	return frequencies / sum(frequencies)
end

# Parameters:
# 	Normal:
# 	- frequencies, integer array
# 	Optional:
# 	- lambda, number OR void
#	- get_target, function
# 	Keyword:
function get_probabilities_shrinkage(frequencies, lambda::Void, get_target = get_uniform_distribution)
	target = get_target(frequencies)
	n = sum(frequencies)
	normalized_frequencies = frequencies / n
	lambda = get_lambda(normalized_frequencies, target, n)
	return apply_shrinkage_formula(normalized_frequencies, target, lambda)
end
function get_probabilities_shrinkage(frequencies, lambda::Number, get_target = get_uniform_distribution)
	target = get_target(frequencies)
	normalized_frequencies = get_normalized_frequencies(frequencies)
	return apply_shrinkage_formula(normalized_frequencies, target, lambda)
end
function apply_shrinkage_formula(normalized_frequencies, target, lambda)
	return lambda * target + (1 - lambda) * normalized_frequencies
end
function get_uniform_distribution(frequencies)
	return 1 / length(frequencies)
end
function get_normalized_frequencies(frequencies)
	return frequencies / sum(frequencies)
end
function get_lambda(normalized_frequencies, target, n)
	if n == 1 || n == 0
		return 1
	end
	# Add better comments about varu and msp
	# Unbiased estimator of variance of u
	varu = normalized_frequencies .* (1 - normalized_frequencies) / (n - 1)
	msp = sum((normalized_frequencies - target).^2) # misspecification ???

	# Estimate shrinkage intensity
	lambda = msp == 0 ? 1 : sum(varu) / msp
	
	# Make lambda be between 0 and 1 inclusive
	return lambda > 1 ? 1 : (lambda < 0 ? 0 : lambda)

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