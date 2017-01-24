# Functions that implement formulae for the information measures

# Information measures are reviewed in:
# Timme, Nicholas; Alford, Wesley; Flecker, Benjamin; Beggs, John M. (2013-07-03).
# "Synergy, redundancy, and multivariate information measures: an experimentalist's perspective".
# Journal of Computational Neuroscience. 36 (2): 119–140.
# http://link.springer.com/article/10.1007%2Fs10827-013-0458-4

# TODO: Add further measures (dual total correlation, delta i)

export apply_entropy_formula, apply_conditional_entropy_formula, apply_mutual_information_formula,
	apply_conditional_mutual_information_formula, apply_interaction_information_formula,
	apply_total_correlation_formula

# Parameters:
# 	- probabilities, array of floats
# 	- base, number
function apply_entropy_formula(probabilities, base)
	probs = probabilities .* log(base, probabilities)
	probs[isnan(probs)] = 0.0
	return -sum(probs)
end

# Parameters:
# 	- joint entropy of all variables, number
#	- entropy of conditioned variables, number
function apply_conditional_entropy_formula(entropy_xy, entropy_y)
	return entropy_xy - entropy_y
end

# Parameters:
# 	- entropy of first variable, number
# 	- entropy of second variable, number
#	- joint entropy of both variables, number
function apply_mutual_information_formula(entropy_x, entropy_y, entropy_xy)
	return entropy_x + entropy_y - entropy_xy
end

# Parameters:
# 	- joint entropy of first and conditioned variables, number
# 	- joint entropy of second and conditioned variables, number
#	- joint entropy of all three variables, number
function apply_conditional_mutual_information_formula(entropy_xz, entropy_yz, entropy_xyz, entropy_z)
	return entropy_xz + entropy_yz - entropy_xyz - entropy_z
end

# Parameters:
# 	- mutual information of first two variables, number
# 	- conditional mutual information of first two variables on the third, number
function apply_interaction_information_formula(conditional_mutual_information, mutual_information)
	return conditional_mutual_information - mutual_information
end

# Parameters:
# 	- entropy of first variable, number
# 	- entropy of second variable, number
# 	- entropy of third variable, number
#	- joint entropy of all three variables, number
function apply_total_correlation_formula(entropy_x, entropy_y, entropy_z, entropy_xyz)
	return entropy_x + entropy_y + entropy_z - entropy_xyz
end
