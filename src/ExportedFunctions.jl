# Functions that the user can call. Should cover all dimensionalities up to 3D:
# discretize
# entropy
# conditional entropy
# mutual information
# conditional mutual information
# total correlation
# dual total correlation
# interaction information
# delta i
# partial information decomposition

# TODO: Credits and citations
# TODO: Add aliases for the modes
# TODO: Add aliases for the estimators
# TODO: Work out whether/how to support values/frequencies/probabilities as user input
# TODO: Add functions for interaction information etc
# TODO: Documentation using Docile and Lexicon

export discretize_values,
	get_probabilities,
	get_entropy,
	get_conditional_entropy,
	get_mutual_information,
	get_conditional_mutual_information,
	get_total_correlation,
	get_dual_total_correlation,
	get_interaction_information,
	get_partial_information_decomposition,
	get_delta_i

# Parameters:
# 	Normal:
#	Varargs:
# 	- values, arrays of floats (could be 1, 2 or 3 - limited by get_frequencies)
# 	Optional:
# 	Keyword:
# 	- mode, string
#	- number_of_bins = 0, number - WARNING: This will get overridden by bayesian_blocks
#	- get_number_of_bins = get_root_n, function - will only be called if number_of_bins is 0
function discretize_values(values...; mode = "uniform_width", number_of_bins = 0, get_number_of_bins = get_root_n)

	if number_of_bins == 0
		number_of_bins = get_number_of_bins(values...)
	end

	return get_frequencies(values..., mode, number_of_bins)
end

# Parameters:
# 	Normal:
# 	- estimator, string
# 	- frequencies, array
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
# 	Keyword:
function get_probabilities(estimator, frequencies, lambda, prior)

	if estimator == "maximum_likelihood" || estimator == "miller_madow"
		probabilities = get_probabilities_maximum_likelihood(frequencies)
	elseif estimator == "shrinkage"
		probabilities = get_probabilities_shrinkage(frequencies, lambda)
	elseif estimator == "dirichlet"
		probabilities = get_probabilities_dirichlet(frequencies, prior)
	end

	return probabilities
end

# Parameters:
# 	Normal:
#	Varargs:
# 	- values, arrays of floats (could be 1, 2 or 3)
# 	Optional:
# 	Keyword:
# 	- estimator, string
# 	- base, number
# 	- mode, string
#	- number_of_bins, number
#	- get_number_of_bins, function
#	- discretized, boolean
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_entropy(values...; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, discretized = false, lambda = nothing, prior = 0)

	# If discretized, values should be one array of frequencies
	# (because 2D frequencies cannot be reconstructed from multiple
	# 1D frequencies), so values... will be ([f1, f2, f3, ...],)
	frequencies = discretized ? values[1] : discretize_values(values..., mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	probabilities = get_probabilities(estimator, frequencies, lambda, prior)
	
	entropy = apply_entropy_formula(probabilities, base)

	if estimator == "miller_madow"
		entropy += (countnz(probabilities) - 1) / (2 * length(probabilities))
	end
	
	return entropy
end

# Parameters:
# 	Normal:
# 	- values_x, array of floats
#	- values_y to condition on
# 	Optional:
# 	Keyword:
# 	- estimator, string
# 	- base, number
# 	- mode, string
#	- number_of_bins, number
#	- get_number_of_bins, function
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_conditional_entropy(values_x, values_y; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 0)
	
	frequencies_xy = discretize_values(values_x, values_y, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_conditional_entropy(frequencies_xy; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
# Parameters:
#	Normal:
#	- xy, array - 2D frequencies or probabilities
#	Keyword:
# 	- estimator, string
# 	- base, number
#	- probabilities = false, boolean - whether xy are probabilities
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_conditional_entropy(xy; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 0)

	probabilities_xy = probabilities ? xy : get_probabilities(estimator, xy, lambda, prior)
	probabilities_y = sum(probabilities_xy, 2)

	entropy_xy = apply_entropy_formula(probabilities_xy, base)
	entropy_y = apply_entropy_formula(probabilities_y, base)

	if estimator == "miller_madow"
		entropy_xy += (countnz(probabilities_xy) - 1) / (2 * length(probabilities_xy))
		entropy_y += (countnz(probabilities_y) - 1) / (2 * length(probabilities_y))
	end

	return apply_conditional_entropy_formula(entropy_xy, entropy_y)
end

# Parameters:
# 	Normal:
# 	- first set of values, array of floats
#	- second set of values
# 	Optional:
# 	Keyword:
# 	- estimator, string
# 	- base, number
# 	- mode, string
#	- number_of_bins, number
#	- get_number_of_bins, function
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_mutual_information(values_x, values_y; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 0)

	frequencies_xy = discretize_values(values_x, values_y, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_mutual_information(frequencies_xy; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
# Parameters:
#	Normal:
#	- xy, array - 2D frequencies or probabilities
#	Keyword:
# 	- estimator, string
# 	- base, number
#	- probabilities = false, boolean - whether xy are probabilities
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_mutual_information(xy; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 0)

	probabilities_xy = probabilities ? xy : get_probabilities(estimator, xy, lambda, prior)
	probabilities_x = sum(probabilities_xy, 1)
	probabilities_y = sum(probabilities_xy, 2)

	entropy_xy = apply_entropy_formula(probabilities_xy, base)
	entropy_x = apply_entropy_formula(probabilities_x, base)
	entropy_y = apply_entropy_formula(probabilities_y, base)

	if estimator == "miller_madow"
		entropy_xy += (countnz(probabilities_xy) - 1) / (2 * length(probabilities_xy))
		entropy_x += (countnz(probabilities_x) - 1) / (2 * length(probabilities_x))
		entropy_y += (countnz(probabilities_y) - 1) / (2 * length(probabilities_y))
	end

	return apply_mutual_information_formula(entropy_x, entropy_y, entropy_xy)
end

# Parameters:
# 	Normal:
# 	- first set of values, array of floats (only 1)
#	- second set of values (only 1)
#	- values to condition on (only 1)
# 	- estimator, string
# 	Optional:
# 	Keyword:
# 	- base, number
# 	- mode, string
#	- number_of_bins, number
#	- get_number_of_bins, function
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_conditional_mutual_information(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 0)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_conditional_mutual_information(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
# Parameters:
#	Normal:
#	- xyz, array - 3D frequencies or probabilities
#	Keyword:
# 	- estimator, string
# 	- base, number
#	- probabilities = false, boolean - whether xyz are probabilities
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_conditional_mutual_information(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 0)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda, prior)
	probabilities_xz = sum(probabilities_xyz, 2)
	probabilities_yz = sum(probabilities_xyz, 1)
	probabilities_z = sum(probabilities_xz, 1) # xz not a typo

	entropy_xyz = apply_entropy_formula(probabilities_xyz, base)
	entropy_xz = apply_entropy_formula(probabilities_xz, base)
	entropy_yz = apply_entropy_formula(probabilities_yz, base)
	entropy_z = apply_entropy_formula(probabilities_z, base)

	if estimator == "miller_madow"
		entropy_xyz += (countnz(probabilities_xyz) - 1) / (2 * length(probabilities_xyz))
		entropy_xz += (countnz(probabilities_xz) - 1) / (2 * length(probabilities_xz))
		entropy_yz += (countnz(probabilities_yz) - 1) / (2 * length(probabilities_yz))
		entropy_z += (countnz(probabilities_z) - 1) / (2 * length(probabilities_z))
	end

	return apply_conditional_mutual_information_formula(entropy_xz, entropy_yz, entropy_xyz, entropy_z)
end

# Parameters:
# 	Normal:
# 	- first set of values, array of floats (only 1)
#	- second set of values (only 1)
#	- values to condition on (only 1)
# 	- estimator, string
# 	Optional:
# 	Keyword:
# 	- base, number
# 	- mode, string
#	- number_of_bins, number
#	- get_number_of_bins, function
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_interaction_information(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 0)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_interaction_information(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
# Parameters:
#	Normal:
#	- xyz, array - 3D frequencies or probabilities
#	Keyword:
# 	- estimator, string
# 	- base, number
#	- probabilities = false, boolean - whether xyz are probabilities
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
function get_interaction_information(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 0)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda, prior)
	probabilities_xy = probabilities_xz = sum(probabilities_xyz, 3)

	conditional_mutual_information = get_conditional_mutual_information(probabilities_xyz, estimator = estimator, base = base, probabilities = true, lambda = lambda, prior = prior)
	mutual_information = get_mutual_information(probabilities_xy, estimator = estimator, base = base, probabilities = true, lambda = lambda, prior = prior)

	return apply_interaction_information_formula(conditional_mutual_information, mutual_information)
end

# Parameters:
# 	Normal:
# 	- first set of values, array of floats (only 1)
#	- second set of values (only 1)
#	- values to condition on (only 1)
# 	- estimator, string
# 	Optional:
# 	Keyword:
# 	- base, number
# 	- mode, string
#	- number_of_bins, number
#	- get_number_of_bins, function
# 	- lambda = nothing, number - only used if estimator is "shrinkage"
# 	- prior = 0, number - only used if estimator is "dirichlet"
#	NB This is only implemented for 3 variables, since:
#	* it is equivalent to mutual information for 2 variables
#	* the package currently doesn't handle more than 3 variables (limited by get_frequencies)
function get_total_correlation(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, discretized = false, lambda = nothing, prior = 0)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_total_correlation(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
function get_total_correlation(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 0)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda, prior)
	probabilities_xy = sum(probabilities_xyz, 3)
	probabilities_x = sum(probabilities_xy, 2)
	probabilities_y = sum(probabilities_xy, 1)
	probabilities_z = sum(sum(probabilities_xyz, 1), 2)

	entropy_xyz = apply_entropy_formula(probabilities_xyz, base)
	entropy_x = apply_entropy_formula(probabilities_x, base)
	entropy_y = apply_entropy_formula(probabilities_y, base)
	entropy_z = apply_entropy_formula(probabilities_z, base)

	return apply_total_correlation_formula(entropy_x, entropy_y, entropy_z, entropy_xyz)
end

function get_partial_information_decomposition()
	function get_redundancy()
		function get_specific_information()
		end
	end
end

function get_dual_total_correlation()
end

function get_delta_i()
end

