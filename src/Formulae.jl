# Functions that apply the formulae for the various information measures

# TODO: use mutliple dispatch so it's possible to call these using the probabilities or already-calculated entropies?
# TODO: generalise to Q formulae and T formulae?
# TODO: DTC and Delta i

export apply_entropy_formula, apply_conditional_entropy_formula, apply_mutual_information_formula,
	apply_conditional_mutual_information_formula, apply_interaction_information_formula,
	apply_total_correlation_formula, apply_dual_total_correlation_formula, apply_delta_i_formula

# Parameters:
# 	Normal:
# 	- probabilities, array of floats
# 	- base, number
# 	Optional:
# 	Keyword:
function apply_entropy_formula(probabilities, base)
	probs = probabilities .* log(base, probabilities)
	probs[isnan(probs)] = 0
	return -sum(probs)
end

# Parameters:
# 	Normal:
# 	- joint entropy of all variables, number
#	- entropy of conditioned variables, number
# 	- base, number
# 	Optional:
# 	Keyword:
function apply_conditional_entropy_formula(entropy_xy, entropy_y)
	return entropy_xy - entropy_y
end

# Parameters:
# 	Normal:
# 	- entropy of first variable, number
# 	- entropy of second variable, number
#	- joint entropy of both variables, number
# 	- base, number
# 	Optional:
# 	Keyword:
function apply_mutual_information_formula(entropy_x, entropy_y, entropy_xy)
	return entropy_x + entropy_y - entropy_xy
end

# Parameters:
# 	Normal:
# 	- joint entropy of first and conditioned variables, number
# 	- joint entropy of second and conditioned variables, number
#	- joint entropy of all three variables, number
# 	- base, number
# 	Optional:
# 	Keyword:
function apply_conditional_mutual_information_formula(entropy_xz, entropy_yz, entropy_xyz, entropy_z)
	return entropy_xz + entropy_yz - entropy_xyz - entropy_z
end

# Parameters:
# 	Normal:
# 	- mutual information of first two variables, number
# 	- conditional mutual information of first two variables on the third, number
# 	- base, number
# 	Optional:
# 	Keyword:
function apply_interaction_information_formula(conditional_mutual_information, mutual_information)
	return conditional_mutual_information - mutual_information
end

# Parameters:
# 	Normal:
# 	- entropy of first variable, number
# 	- entropy of second variable, number
# 	- entropy of third variable, number
#	- joint entropy of all three variables, number
# 	- base, number
# 	Optional:
# 	Keyword:
function apply_total_correlation_formula(entropy_x, entropy_y, entropy_z, entropy_xyz)
	return entropy_x + entropy_y + entropy_z - entropy_xyz
end
