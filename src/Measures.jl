# Functions for estimating information measures for continuous or discrete datasets.
# These are intended to maximise ease of use for setting different options. For
# maximum speed, it may be better to use the functions from Formulae.jl directly.

export discretize_values,
	get_probabilities,
	get_entropy,
	get_conditional_entropy,
	get_mutual_information,
	get_conditional_mutual_information,
	get_total_correlation,
	get_interaction_information,
	get_partial_information_decomposition,
	get_redundancy,
	get_cross_entropy

"""
    discretize_values(values_x; <keyword arguments>)

Assign continuous measurements to discrete bins.

# Arguments
* `values_x`: the data values.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
"""
function discretize_values(values_x...; mode = "uniform_width", number_of_bins = 0, get_number_of_bins = get_root_n)

	if number_of_bins == 0
		number_of_bins = get_number_of_bins(values_x...)
	end

	return get_frequencies(mode, number_of_bins, values_x...)
end

"""
    get_probabilities(estimator, frequencies; <keyword arguments>)

Estimate probabilities from a set of discrete values.

# Arguments:
* `estimator`: the entropy estimator.
* `frequencies`: the bin frequencies for the discretized data values.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_probabilities(estimator, frequencies; lambda = nothing, prior = 1)

	if estimator == "maximum_likelihood" || estimator == "miller_madow"
		probabilities = get_probabilities_maximum_likelihood(frequencies)
	elseif estimator == "shrinkage"
		probabilities = get_probabilities_shrinkage(frequencies, lambda)
	elseif estimator == "dirichlet"
		probabilities = get_probabilities_dirichlet(frequencies, prior)
	end

	return probabilities
end

"""
    get_entropy(values_x; <keyword arguments>)

Estimate entropy (or joint entropy, if more than one set of data values is given).

# Arguments:
* `values_x`: the data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `discretized=false`: whether the data values are already discretized.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_entropy(values_x...; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, discretized = false, lambda = nothing, prior = 1)

	# If discretized, values_x should be one array of frequencies
	# (because higher dimensional frequencies cannot be reconstructed
	# from multiple lower dimensional frequencies), so values_x... will
	# be ([f1, f2, f3, ...],)
	frequencies = discretized ?
		values_x[1] :
		discretize_values(values_x..., mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	probabilities = get_probabilities(estimator, frequencies, lambda = lambda, prior = prior)

	entropy = apply_entropy_formula(probabilities, base)

	if estimator == "miller_madow"
		entropy += (count(!iszero, probabilities) - 1) / (2 * length(values_x[1]))
	end

	# Make sure we don't return -0.0
	return entropy + 0.0
end

"""
    get_conditional_entropy(values_x, values_y; <keyword arguments>)

Estimate conditional entropy of one set of values conditioned on another set of values.

# Arguments:
* `values_x`: the data values.
* `values_y`: the data values to be conditioned on.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_conditional_entropy(values_x, values_y; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 1)

	frequencies_xy = discretize_values(values_x, values_y, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_conditional_entropy(frequencies_xy; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
"""
    get_conditional_entropy(xy; <keyword arguments>)

# Arguments:
* `xy`: the 2D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xy` is probabilities. If `false`, `xy` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_conditional_entropy(xy; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1)

	probabilities_xy = probabilities ? xy : get_probabilities(estimator, xy, lambda = lambda, prior = prior)
	probabilities_y = sum(probabilities_xy, dims = 2)

	entropy_xy = apply_entropy_formula(probabilities_xy, base)
	entropy_y = apply_entropy_formula(probabilities_y, base)

	# TODO: either deprecate or add warning about inappropriate use
	if estimator == "miller_madow"
		entropy_xy += (count(!iszero, probabilities_xy) - 1) / (2 * length(probabilities_xy))
		entropy_y += (count(!iszero, probabilities_y) - 1) / (2 * length(probabilities_y))
	end

	return apply_conditional_entropy_formula(entropy_xy, entropy_y)
end

"""
    get_mutual_information(values_x, values_y; <keyword arguments>)

Estimate the mutual information between two sets of values.

# Arguments:
* `values_x`: the data values.
* `values_y`: a second set of data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_mutual_information(values_x, values_y; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 1)

	frequencies_xy = discretize_values(values_x, values_y, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_mutual_information(frequencies_xy; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
"""
    get_mutual_information(xy; <keyword arguments>)

# Arguments:
* `xy`: the 2D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xy` is probabilities. If `false`, `xy` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_mutual_information(xy; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1)

	probabilities_xy = probabilities ? xy : get_probabilities(estimator, xy, lambda = lambda, prior = prior)
	probabilities_x = sum(probabilities_xy, dims = 1)
	probabilities_y = sum(probabilities_xy, dims = 2)

	# TODO: either deprecate or add warning about inappropriate use
	if estimator == "miller_madow"
		entropy_xy = apply_entropy_formula(probabilities_xy, base)
		entropy_x = apply_entropy_formula(probabilities_x, base)
		entropy_y = apply_entropy_formula(probabilities_y, base)
		entropy_xy += (count(!iszero, probabilities_xy) - 1) / (2 * length(probabilities_xy))
		entropy_x += (count(!iszero, probabilities_x) - 1) / (2 * length(probabilities_x))
		entropy_y += (count(!iszero, probabilities_y) - 1) / (2 * length(probabilities_y))
		return apply_mutual_information_formula(entropy_x, entropy_y, entropy_xy)
	end

	return apply_mutual_information_formula(probabilities_xy, probabilities_x, probabilities_y, base)
end

"""
    get_conditional_mutual_information(values_x, values_y, values_z; <keyword arguments>)

Estimate the conditional mutual information between two sets of values, conditioned on a third.

# Arguments:
* `values_x`: the data values.
* `values_y`: a second set of data values.
* `values_z`: the data values to be conditioned on.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_conditional_mutual_information(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 1)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_conditional_mutual_information(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
"""
    get_conditional_mutual_information(xyz; <keyword arguments>)

# Arguments:
* `xyz`: the 3D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xyz` is probabilities. If `false`, `xyz` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_conditional_mutual_information(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda = lambda, prior = prior)
	probabilities_xz = sum(probabilities_xyz, dims = 2)
	probabilities_yz = sum(probabilities_xyz, dims = 1)
	probabilities_z = sum(probabilities_xz, dims = 1) # xz not a typo

	entropy_xyz = apply_entropy_formula(probabilities_xyz, base)
	entropy_xz = apply_entropy_formula(probabilities_xz, base)
	entropy_yz = apply_entropy_formula(probabilities_yz, base)
	entropy_z = apply_entropy_formula(probabilities_z, base)

	# TODO: either deprecate or add warning about inappropriate use
	if estimator == "miller_madow"
		entropy_xyz += (count(!iszero, probabilities_xyz) - 1) / (2 * length(probabilities_xyz))
		entropy_xz += (count(!iszero, probabilities_xz) - 1) / (2 * length(probabilities_xz))
		entropy_yz += (count(!iszero, probabilities_yz) - 1) / (2 * length(probabilities_yz))
		entropy_z += (count(!iszero, probabilities_z) - 1) / (2 * length(probabilities_z))
	end

	return apply_conditional_mutual_information_formula(entropy_xz, entropy_yz, entropy_xyz, entropy_z)
end

"""
    get_interaction_information(values_x, values_y, values_z; <keyword arguments>)

Estimate the interaction information between three sets of values.

# Arguments:
* `values_x`: the data values.
* `values_y`: a second set of data values.
* `values_z`: a third set of data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_interaction_information(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width",
	number_of_bins = 0, get_number_of_bins = get_root_n, lambda = nothing, prior = 1)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_interaction_information(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
"""
    get_interaction_information(xyz; <keyword arguments>)

# Arguments:
* `xyz`: the 3D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xyz` is probabilities. If `false`, `xyz` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_interaction_information(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda = lambda, prior = prior)
	probabilities_xy = probabilities_xz = sum(probabilities_xyz, dims = 3)

	conditional_mutual_information = get_conditional_mutual_information(probabilities_xyz, estimator = estimator, base = base, probabilities = true, lambda = lambda, prior = prior)
	mutual_information = get_mutual_information(probabilities_xy, estimator = estimator, base = base, probabilities = true, lambda = lambda, prior = prior)

	return apply_interaction_information_formula(conditional_mutual_information, mutual_information)
end

"""
    get_total_correlation(values_x, values_y, values_z; <keyword arguments>)

Estimate the total correlation between three sets of values.

# Arguments:
* `values_x`: the data values.
* `values_y`: a second set of data values.
* `values_z`: a third set of data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_total_correlation(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, lambda = nothing, prior = 1)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_total_correlation(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior)
end
"""
    get_total_correlation(xyz; <keyword arguments>)

# Arguments:
* `xyz`: the 3D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xyz` is probabilities. If `false`, `xyz` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_total_correlation(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda = lambda, prior = prior)
	probabilities_xy = sum(probabilities_xyz, dims = 3)
	probabilities_x = sum(probabilities_xy, dims = 2)
	probabilities_y = sum(probabilities_xy, dims = 1)
	probabilities_z = sum(sum(probabilities_xyz, dims = 1), dims = 2)

	entropy_xyz = apply_entropy_formula(probabilities_xyz, base)
	entropy_x = apply_entropy_formula(probabilities_x, base)
	entropy_y = apply_entropy_formula(probabilities_y, base)
	entropy_z = apply_entropy_formula(probabilities_z, base)

	return apply_total_correlation_formula(entropy_x, entropy_y, entropy_z, entropy_xyz)
end

"""
    get_partial_information_decomposition(values_x, values_y, values_z; <keyword arguments>)

Estimate the partial information decomposition between three sets of values.

Performance can be improved by setting `all_orientations`, `include_unique` and `include_synergy` according to the use case.

# Arguments:
* `values_x`: the data values.
* `values_y`: a second set of data values.
* `values_z`: a third set of data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equqivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
* `all_orientations=false`: whether to use each set of values as the target. If `false`, only `values_z` is the target.
* `include_unique=true`: whether to get the unique information.
* `include_synergy=true`: whether to get the synergy.
"""
function get_partial_information_decomposition(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, lambda = nothing, prior = 1, all_orientations = false, include_unique = true, include_synergy = true)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_partial_information_decomposition(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior, all_orientations = all_orientations,
		include_unique = include_unique, include_synergy = include_synergy)
end
"""
    get_partial_information_decomposition(xyz; <keyword arguments>)

# Arguments:
* `xyz`: the 3D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xyz` is probabilities. If `false`, `xyz` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
* `all_orientations=false`: whether to use each set of values as the target. If `false`, only `values_z` is the target.
* `include_unique=true`: whether to get the unique information.
* `include_synergy=true`: whether to get the synergy.
"""
function get_partial_information_decomposition(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1, all_orientations = false,
	include_unique = true, include_synergy = true)

	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda = lambda, prior = prior)

	probabilities_xz = sum(probabilities_xyz, dims = 2)
	probabilities_yz = sum(probabilities_xyz, dims = 1)
	probabilities_xy = sum(probabilities_xyz, dims = 3)
	probabilities_x = sum(probabilities_xz, dims = 3)
	probabilities_y = sum(probabilities_yz, dims = 3)
	probabilities_z = sum(probabilities_xz, dims = 1)

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
		redundancy_xy_z = apply_redundancy_formula(probabilities_xz, probabilities_yz, probabilities_x, probabilities_y, probabilities_z, (1, 2), base)
		redundancy_xz_y = apply_redundancy_formula(probabilities_xy, probabilities_yz, probabilities_x, probabilities_z, probabilities_y, (1, 3), base)
		redundancy_yz_x = apply_redundancy_formula(probabilities_xy, probabilities_xz, probabilities_y, probabilities_z, probabilities_x, (2, 3), base)

		pid["z"] = Dict(
			"redundancy" => redundancy_xy_z
		)
		pid["y"] = Dict(
			"redundancy" => redundancy_xz_y
		)
		pid["x"] = Dict(
			"redundancy" => redundancy_yz_x
		)

		if include_unique
			pid["z"]["unique_1"] = mutual_information_xz - redundancy_xy_z
			pid["z"]["unique_2"] = mutual_information_yz - redundancy_xy_z
			pid["y"]["unique_1"] = mutual_information_xy - redundancy_xz_y
			pid["y"]["unique_2"] = mutual_information_yz - redundancy_xz_y
			pid["x"]["unique_1"] = mutual_information_xy - redundancy_yz_x
			pid["x"]["unique_2"] = mutual_information_xz - redundancy_yz_x
		end

		if include_synergy
			pid["z"]["synergy"] = interaction_information + redundancy_xy_z
			pid["y"]["synergy"] = interaction_information + redundancy_xz_y
			pid["x"]["synergy"] = interaction_information + redundancy_yz_x
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

		redundancy = apply_redundancy_formula(probabilities_xz, probabilities_yz, probabilities_x, probabilities_y, probabilities_z, (1, 2), base)
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

"""
    get_redundancy(values_x, values_y, values_z; <keyword arguments>)

Estimate the redundancy between three sets of values.

# Arguments:
* `values_x`: the data values.
* `values_y`: a second set of data values.
* `values_z`: a third set of data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equqivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
* `all_orientations=false`: whether to use each set of values as the target. If `false`, only `values_z` is the target.
"""
function get_redundancy(values_x, values_y, values_z; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0,
	get_number_of_bins = get_root_n, lambda = nothing, prior = 1, all_orientations = false)

	frequencies_xyz = discretize_values(values_x, values_y, values_z, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

	return get_redundancy(frequencies_xyz; estimator = estimator, base = base, lambda = lambda, prior = prior, all_orientations = all_orientations)
end
"""
    get_redundancy(xyz; <keyword arguments>)

# Arguments:
* `xyz`: the 3D frequencies or probabilities.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `probabilities=false`: whether `xyz` is probabilities. If `false`, `xyz` is frequencies.
* `lambda=nothing`: the shrinkage instensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
* `all_orientations=false`: whether to use each set of values as the target. If `false`, only `values_z` is the target.
"""
function get_redundancy(xyz; estimator = "maximum_likelihood", base = 2, probabilities = false, lambda = nothing, prior = 1, all_orientations = false)
	probabilities_xyz = probabilities ? xyz : get_probabilities(estimator, xyz, lambda = lambda, prior = prior)

	probabilities_xz = sum(probabilities_xyz, dims = 2)
	probabilities_yz = sum(probabilities_xyz, dims = 1)
	probabilities_xy = sum(probabilities_xyz, dims = 3)
	probabilities_x = sum(probabilities_xz, dims = 3)
	probabilities_y = sum(probabilities_yz, dims = 3)
	probabilities_z = sum(probabilities_xz, dims = 1)

	if all_orientations
		redundancy_xy_z = apply_redundancy_formula(probabilities_xz, probabilities_yz, probabilities_x, probabilities_y, probabilities_z, (1, 2), base)
		redundancy_xz_y = apply_redundancy_formula(probabilities_xy, probabilities_yz, probabilities_x, probabilities_z, probabilities_y, (1, 3), base)
		redundancy_yz_x = apply_redundancy_formula(probabilities_xy, probabilities_xz, probabilities_y, probabilities_z, probabilities_x, (2, 3), base)

		return (redundancy_xy_z, redundancy_xz_y, redundancy_yz_x)
	else
		return apply_redundancy_formula(probabilities_xz, probabilities_yz, probabilities_x, probabilities_y, probabilities_z, (1, 2), base)
	end
end


"""
    get_cross_entropy(values_x, values_y; <keyword arguments>)
Estimate the cross-entropy between two sets of values.
# Arguments:
* `values_x`: the first set of data values.
* `values_y`: the second set of data values.
* `estimator="maximum_likelihood"`: the entropy estimator.
* `base=2`: the base of the logarithm, equivalent to the units of information.
* `mode="uniform_width"`: the discretization algorithm.
* `number_of_bins=0`: the number of bins (will be overridden if `mode` is `"bayesian_blocks"`).
* `get_number_of_bins=get_root_n`: the method for calculating the number of bins (only called if `number_of_bins` is `0`).
* `discretized=false`: whether the data values are already discretized.
* `lambda=nothing`: the shrinkage intensity, only used if `estimator` is `"shrinkage"`.
* `prior=1`: the Dirichlet prior, only used if `estimator` is `"dirichlet"`.
"""
function get_cross_entropy(values_x, values_y; estimator = "maximum_likelihood", base = 2, mode = "uniform_width", number_of_bins = 0, get_number_of_bins = get_root_n, discretized = false, lambda = nothing, prior = 1)
    frequency_x = discretized ? values_x : discretize_values(values_x, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)
    frequency_y = discretized ? values_y : discretize_values(values_y, mode = mode, number_of_bins = number_of_bins, get_number_of_bins = get_number_of_bins)

    probability_x = get_probabilities(estimator, frequency_x, lambda = lambda, prior = prior)
    probability_y = get_probabilities(estimator, frequency_y, lambda = lambda, prior = prior)

    cross_entropy = apply_cross_entropy_formula(probability_x, probability_y, base)

    if estimator == "miller_madow"
        println("WARNING: Miller-Madow correction not implemented for the cross-entropy. ")
        println("The calculation was performed without applying the correction.")
    end

    return cross_entropy
end
