# Functions that implement formulae for the information measures

# Information measures are reviewed in:
# Timme, Nicholas; Alford, Wesley; Flecker, Benjamin; Beggs, John M. (2013-07-03).
# "Synergy, redundancy, and multivariate information measures: an experimentalist's perspective".
# Journal of Computational Neuroscience. 36 (2): 119–140.
# http://link.springer.com/article/10.1007%2Fs10827-013-0458-4

# TODO: Add further measures (dual total correlation, delta i)

export apply_entropy_formula, apply_conditional_entropy_formula, apply_mutual_information_formula,
	apply_conditional_mutual_information_formula, apply_interaction_information_formula,
	apply_total_correlation_formula, apply_specific_information_formula, apply_redundancy_formula,
	apply_unique_information_formula, apply_synergy_formula, apply_cross_entropy_formula

function remove_non_finite(x)
	return isfinite(x) ? x : zero(x)
end

function remove_negative(x)
	return x < 0 ? zero(x) : x
end

# Parameters:
# 	- probabilities, array of floats
# 	- base, number
function apply_entropy_formula(p::AbstractArray{T}, base::R) where {T <: AbstractFloat, R <: Real}
	return -sum(remove_non_finite.(p .* log.(base, p)))
end

# Parameters:
# 	- joint entropy of all variables, number
#	- entropy of conditioned variables, number
function apply_conditional_entropy_formula(entropy_xy, entropy_y)
	return entropy_xy - entropy_y
end

# Parameters:
# 	- joint probabilities, array of floats
# 	- probabilities (first variable), array of floats
# 	- probabilities (second variable), array of floats
# 	- base, number
function apply_mutual_information_formula(p_xy::AbstractArray{T}, p_x::AbstractArray{T}, p_y::AbstractArray{T}, base::R) where {T <: AbstractFloat, R <: Real}
	return sum(remove_non_finite.(p_xy .* log.(base, p_xy ./ (p_x .* p_y))))
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

# Parameters:
# 	- joint probabilities, array of floats
# 	- probabilities (source), array of floats
# 	- probabilities (target), array of floats
# 	- dimension along which to sum, integer
# 	- base, number
function apply_specific_information_formula(p_xz, p_x, p_z, dim_sum, base)
	return vec(sum(remove_non_finite.((p_xz ./ p_z) .* log.(base, p_xz ./ (p_x .* p_z))), dims = dim_sum))
end

# Parameters:
# 	- joint probabilities (source 1 and target), array of floats
# 	- joint probabilities (source 2 and target), array of floats
# 	- probabilities (source 1), array of floats
# 	- probabilities (source 2), array of floats
# 	- probabilities (target), array of floats
# 	- dimensions along which to sum, tuple of integers
# 	- base, number
function apply_redundancy_formula(p_xz::AbstractArray{T}, p_yz::AbstractArray{T}, p_x::AbstractArray{T}, p_y::AbstractArray{T}, p_z::AbstractArray{T}, dim_sum::Tuple{I, I}, base::R) where {T <: AbstractFloat, R <: Real, I <: Integer}
	minimum_specific_information = min.(
		apply_specific_information_formula(p_xz, p_x, p_z, dim_sum[1], base),
		apply_specific_information_formula(p_yz, p_y, p_z, dim_sum[2], base)
	)
	return sum(vec(p_z) .* vec(minimum_specific_information))
end
# Parameters:
# 	- probabilities (target), array of floats
# 	- specific information of source 1 and target, array of floats
# 	- specific information of source 2 and target, array of floats
# 	- base, number
function apply_redundancy_formula(p_z::AbstractArray{T}, specific_information_1::AbstractArray{T}, specific_information_2::AbstractArray{T}, base::R) where {T <: AbstractFloat, R <: Real}
	minimum_specific_information = min.(specific_information_1, specific_information_2)
	return sum(vec(p_z) .* vec(minimum_specific_information))
end

# Parameters:
# 	- mutual information of source 1 and target, number
# 	- redundancy of both sources and target, number
function apply_unique_information_formula(mutual_information, redundancy)
	# Rounding errors may lead to slightly negative results
	return remove_negative.(mutual_information - redundancy)
end

# Parameters:
# 	- interaction information of both sources and target, number
# 	- redundancy of both sources and target, number
function apply_synergy_formula(interaction_information, redundancy)
	return interaction_information + redundancy
end

# Parameters:
#	- probabilities of first variable
#	- probabilities of second variable
#	- base of exponent
function apply_cross_entropy_formula(p_x::AbstractArray{T}, p_y::AbstractArray{T}, base::R) where {T <: AbstractFloat, R <: Real}
    return -sum(remove_non_finite.(p_x .* log.(base, p_y)))
end
