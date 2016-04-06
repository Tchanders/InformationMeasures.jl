# Functions to aid discretization of raw values

using Discretizers
using PyCall

@pyimport astroML.density_estimation as de

export get_bin_ids!, get_frequencies_from_bin_ids

# Parameters:
# 	Normal:
#	Varargs:
# 	- values, arrays of floats (could be 1, 2 or 3)
# 	Optional:
# 	Keyword:
function get_root_n(values...)
	return round(Int, sqrt(size(values[1])[1]))
end

# Parameters:
# 	Normal:
# 	- values, arrays of floats (could be 1, 2 or 3)
# 	- mode, string
#	- number_of_bins, number
# 	Optional:
# 	Keyword:
function get_frequencies(mode, number_of_bins, values_x; discretized = false)
	bin_ids = Array(Int, size(values_x))
	number_of_bins = get_bin_ids!(values_x, mode, number_of_bins, bin_ids)
	frequencies = zeros(Int, (number_of_bins, 1))
	for id in bin_ids
		frequencies[id] += 1
	end
	return frequencies
end
function get_frequencies(mode, number_of_bins, values_x, values_y)
	# n is the number of data points
	n = size(values_x)[1]
	bin_ids_x = Array(Int, size(values_x))
	bin_ids_y = Array(Int, size(values_y))
	number_of_bins_x = get_bin_ids!(values_x, mode, number_of_bins, bin_ids_x)
	number_of_bins_y = get_bin_ids!(values_y, mode, number_of_bins, bin_ids_y)
	frequencies = zeros(Int, (number_of_bins_x, number_of_bins_y))
	for i in 1:n
		frequencies[bin_ids_x[i], bin_ids_y[i]] += 1
	end
	return frequencies
end
function get_frequencies(mode, number_of_bins, values_x, values_y, values_z)
	# n is the number of data points
	n = size(values_x)[1]
	bin_ids_x = Array(Int, size(values_x))
	bin_ids_y = Array(Int, size(values_y))
	bin_ids_z = Array(Int, size(values_z))
	number_of_bins_x = get_bin_ids!(values_x, mode, number_of_bins, bin_ids_x)
	number_of_bins_y = get_bin_ids!(values_y, mode, number_of_bins, bin_ids_y)
	number_of_bins_z = get_bin_ids!(values_z, mode, number_of_bins, bin_ids_z)
	frequencies = zeros(Int, (number_of_bins_x, number_of_bins_y, number_of_bins_z))
	for i in 1:n
		frequencies[bin_ids_x[i], bin_ids_y[i], bin_ids_z[i]] += 1
	end
	return frequencies
end
# TODO: sort out the ordering of the parameters
function get_frequencies(mode, number_of_bins, values...)
	# n is the number of data points
	n = size(values[1])[1]
	# d is the number of dimensions
	d = length(values)
	bin_ids = Array(Array{Int}, d)
	all_number_of_bins = Array(Int, d)
	for i in 1:d
		current_values = values[i]
		bin_ids[i] = Array(Int, size(current_values))
		all_number_of_bins[i] = get_bin_ids!(current_values, mode, number_of_bins, bin_ids[i])
	end
	println(all_number_of_bins)
	frequencies = zeros(Int, tuple(all_number_of_bins...))
	for i in 1:n
		frequencies[ntuple(x -> bin_ids[x][i], d)...] += 1
	end
	return frequencies
end

# Added this function to speed up PIDs for large networks
function get_frequencies_from_bin_ids(bin_ids_x, bin_ids_y, number_of_bins_x, number_of_bins_y)
	n = size(bin_ids_x)[1]
	frequencies = zeros(Int, (number_of_bins_x, number_of_bins_y))
	for i in 1:n
		frequencies[bin_ids_x[i], bin_ids_y[i]] += 1
	end
	return frequencies
end
function get_frequencies_from_bin_ids(bin_ids_x, bin_ids_y, bin_ids_z, number_of_bins_x, number_of_bins_y, number_of_bins_z)
	n = size(bin_ids_x)[1]
	frequencies = zeros(Int, (number_of_bins_x, number_of_bins_y, number_of_bins_z))
	for i in 1:n
		frequencies[bin_ids_x[i], bin_ids_y[i], bin_ids_z[i]] += 1
	end
	return frequencies
end

# Parameters:
# 	Normal:
# 	- values, arrays of floats (coule be 1, 2 or 3)
# 	- mode, number
#	- number_of_bins, integer
#	- bin_ids, 1-dimensional array of bin ids for each value
# 	Optional:
# 	Keyword:
# function get_bin_ids(values, mode, number_of_bins)
function get_bin_ids!(values, mode, number_of_bins, bin_ids)
	min, max = extrema(values)
	if min == max
		# If values are all the same, assign them all to bin 1
		number_of_bins = 1
		bin_ids[1:end] = 1
	elseif mode == "uniform_width"
		bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values)), values)
	elseif mode == "binarize"
		number_of_bins = 2
		bin_ids[1:end] = map(v -> v == 0 ? 1 : 2, values)
	elseif mode == "uniform_count"
		try
			bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformCount(number_of_bins), reshape(values, length(values)))), values)
		catch
			bin_ids += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values)), values)
			println("Uniform count failed, fell back to uniform width")
		end
	elseif mode == "bayesian_blocks"
		try
			# The commented line is for a future Julia implementation of Bayesian blocks
			# edges = sort(unique(bayesianBlocks(reshape(values, length(values)))))
			edges = de.bayesian_blocks(reshape(values, length(values)))
			bin_ids[1:end] = encode(LinearDiscretizer(edges), values)
			number_of_bins = length(edges) - 1
		catch
			bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values)), values)
			println("Bayesian blocks failed, fell back to uniform width")
		end
	else
		println("Mode doesn't exist, fell back to uniform width")
		bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values)), values)
	end

	return number_of_bins
end