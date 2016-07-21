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
# 	- values_x, arrays of floats (could be 1, 2 or 3 - limited by get_frequencies)
# 	Optional:
# 	Keyword:
# 	- mode, string
#	- number_of_bins = 0, number - WARNING: This will get overridden by bayesian_blocks
#	- get_number_of_bins = get_root_n, function - will only be called if number_of_bins is 0
function discretize_values(values_x...; mode = "uniform_width", number_of_bins = 0, get_number_of_bins = get_root_n)

	if number_of_bins == 0
		number_of_bins = get_number_of_bins(values_x...)
	end

	if length(values_x) > 3
		return get_frequencies(mode, number_of_bins, values_x...)
	end

	return get_frequencies(mode, number_of_bins, values_x...)
end

# TODO: Change order of estimator and frequencies?
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
# 	- values_x, arrays of floats (could be 1, 2 or 3)
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
function get_entropy(values_x...; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, discretized = false, lambda = nothing, prior = 0)

	# If discretized, values_x should be one array of frequencies
	# (because higher dimensional frequencies cannot be reconstructed
	# from multiple lower dimensional frequencies), so values_x... will
	# be ([f1, f2, f3, ...],)
	frequencies = discretized ?
		values_x[1] :
		discretize_values(values_x..., mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	probabilities = get_probabilities(estimator, frequencies, lambda, prior)
	
	entropy = apply_entropy_formula(probabilities, base)

	if estimator == "miller_madow"
		entropy += (countnz(probabilities) - 1) / (2 * length(values_x[1]))
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

	# TODO: correct this (see get_entropy)
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

	# TODO: correct this (see get_entropy)
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

	# TODO: correct this (see get_entropy)
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

# TODO: add documentation
# TODO: Add Miller-Madow correction
# TODO: Refactor so this always returns the same type
# - all_orientations, boolean - whether or not to calculate pid with each of the three variables as targets
function get_partial_information_decomposition(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, discretized = false, lambda = nothing, prior = 0, all_orientations = false, include_unique = true, include_synergy = true)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_partial_information_decomposition(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior, all_orientations = all_orientations,
		include_unique = include_unique, include_synergy = include_synergy)
end

function get_partial_information_decomposition(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 0, all_orientations = false,
	include_unique = true, include_synergy = true)

	function get_redundancy(probabilities_target, probabilities_source_1, probabilities_source_2, probabilities_source_1_target, probabilities_source_2_target, target_dimension)

		function get_specific_information(probabilities_source_target_i, probabilities_source, probabilities_target_i)
			specific_information = (probabilities_source_target_i / probabilities_target_i) .* log(base, probabilities_source_target_i ./ (probabilities_source * probabilities_target_i))
			specific_information[!isfinite(specific_information)] = 0
			return sum(specific_information)
		end

		number_of_target_states = length(probabilities_target)
		minimum_specific_information = zeros(number_of_target_states)

		for (i, probability_target) in enumerate(probabilities_target)
			if target_dimension == 3
				specific_information_source_1 = get_specific_information(probabilities_source_1_target[:, :, i], probabilities_source_1, probability_target)
				specific_information_source_2 = get_specific_information(probabilities_source_2_target[:, :, i], probabilities_source_2, probability_target)
			elseif target_dimension == 2
				specific_information_source_1 = get_specific_information(probabilities_source_1_target[:, i, :], probabilities_source_1, probability_target)
				specific_information_source_2 = get_specific_information(probabilities_source_2_target[:, i, :], probabilities_source_2, probability_target)
			else
				specific_information_source_1 = get_specific_information(probabilities_source_1_target[i, :, :], probabilities_source_1, probability_target)
				specific_information_source_2 = get_specific_information(probabilities_source_2_target[i, :, :], probabilities_source_2, probability_target)
			end
			minimum_specific_information[i] += min(specific_information_source_1, specific_information_source_2)
		end

		return sum(collect(probabilities_target) .* minimum_specific_information)
	end

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda, prior)

	probabilities_xz = sum(probabilities_xyz, 2)
	probabilities_yz = sum(probabilities_xyz, 1)
	probabilities_xy = sum(probabilities_xyz, 3)
	probabilities_x = sum(probabilities_xz, 3)
	probabilities_y = sum(probabilities_yz, 3)
	probabilities_z = sum(probabilities_xz, 1)

	# Needed for unique information and synergy
	if include_unique || include_synergy
		entropy_x = apply_entropy_formula(probabilities_x, base)
		entropy_y = apply_entropy_formula(probabilities_y, base)
		entropy_z = apply_entropy_formula(probabilities_z, base)
		entropy_xy = apply_entropy_formula(probabilities_xy, base)
		entropy_xz = apply_entropy_formula(probabilities_xz, base)
		entropy_yz = apply_entropy_formula(probabilities_yz, base)
		mutual_information_xy = apply_mutual_information_formula(entropy_x, entropy_y, entropy_xy)
	end

	# Needed for unique information only
	if include_unique
		mutual_information_xz = apply_mutual_information_formula(entropy_x, entropy_z, entropy_xz)
		mutual_information_yz = apply_mutual_information_formula(entropy_y, entropy_z, entropy_yz)
	end

	# Needed for synergy only
	if include_synergy
		entropy_xyz = apply_entropy_formula(probabilities_xyz, base)
		interaction_information = apply_interaction_information_formula(
			apply_conditional_mutual_information_formula(entropy_xz, entropy_yz, entropy_xyz, entropy_z), mutual_information_xy
		)
	end

	pid = Dict()

	if all_orientations
		redundnacy_xy_z = get_redundancy(probabilities_z, probabilities_x, probabilities_y, probabilities_xz, probabilities_yz, 3)
		redundnacy_xz_y = get_redundancy(probabilities_y, probabilities_x, probabilities_z, probabilities_xy, probabilities_yz, 2)
		redundnacy_yz_x = get_redundancy(probabilities_x, probabilities_y, probabilities_z, probabilities_xy, probabilities_xz, 1)

		pid["z"] = Dict(
			"redundancy" => redundnacy_xy_z
		)
		pid["y"] = Dict(
			"redundancy" => redundnacy_xz_y
		)
		pid["x"] = Dict(
			"redundancy" => redundnacy_yz_x
		)

		if include_unique
			pid["z"]["unique_1"] = mutual_information_xz - redundnacy_xy_z
			pid["z"]["unique_2"] = mutual_information_yz - redundnacy_xy_z
			pid["y"]["unique_1"] = mutual_information_xy - redundnacy_xz_y
			pid["y"]["unique_2"] = mutual_information_yz - redundnacy_xz_y
			pid["x"]["unique_1"] = mutual_information_xy - redundnacy_yz_x
			pid["x"]["unique_2"] = mutual_information_xz - redundnacy_yz_x
		end

		if include_synergy
			pid["z"]["synergy"] = interaction_information + redundnacy_xy_z
			pid["y"]["synergy"] = interaction_information + redundnacy_xz_y
			pid["x"]["synergy"] = interaction_information + redundnacy_yz_x
		end

		# Rounding errors may lead to slightly negative results
		for orientation in ["x", "y", "z"]
			for measure in ["redundancy", "unique_1", "unique_2", "synergy"]
				if haskey(pid[orientation], measure)
					pid[orientation][measure] = pid[orientation][measure] < 0 ? 0 : pid[orientation][measure]
				end
			end
		end

	else
		redundancy = get_redundancy(probabilities_z, probabilities_x, probabilities_y, probabilities_xz, probabilities_yz, 3)
		pid["redundancy"] = redundancy

		if include_unique
			unique_x = mutual_information_xz - redundancy
			unique_y = mutual_information_yz - redundancy
			# Rounding errors may lead to slightly negative results
			pid["unique_1"] = unique_x < 0 ? 0 : unique_x
			pid["unique_2"] = unique_y < 0 ? 0 : unique_y
		end

		if include_synergy
			synergy = interaction_information + redundancy
			# Rounding errors may lead to slightly negative results
			pid["synergy"] = synergy < 0 ? 0 : synergy
		end
	end

	return pid
end

function get_dual_total_correlation()
end

function get_delta_i()
end

