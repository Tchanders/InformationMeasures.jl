# Functions to ia ddiscretization of raw values

using Discretizers

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
function get_frequencies(values_x, mode, number_of_bins)
	bin_ids, number_of_bins = get_bin_ids(values_x, mode, number_of_bins)
	frequencies = zeros(Int, (number_of_bins, 1))
	for id in bin_ids
		frequencies[id] += 1
	end
	return frequencies
end
function get_frequencies(values_x, values_y, mode, number_of_bins)
	# n is the number of data points
	n = size(values_x)[1]
	bin_ids_x, number_of_bins_x = get_bin_ids(values_x, mode, number_of_bins)
	bin_ids_y, number_of_bins_y = get_bin_ids(values_y, mode, number_of_bins)
	frequencies = zeros(Int, (number_of_bins_x, number_of_bins_y))
	for i in 1:n
		frequencies[bin_ids_x[i], bin_ids_y[i]] += 1
	end
	return frequencies
end
function get_frequencies(values_x, values_y, values_z, mode, number_of_bins)
	# n is the number of data points
	n = size(values_x)[1]
	bin_ids_x, number_of_bins_x = get_bin_ids(values_x, mode, number_of_bins)
	bin_ids_y, number_of_bins_y = get_bin_ids(values_y, mode, number_of_bins)
	bin_ids_z, number_of_bins_z = get_bin_ids(values_z, mode, number_of_bins)
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
# 	Optional:
# 	Keyword:
function get_bin_ids(values, mode, number_of_bins)
	# bin_ids is a 1-dimensional array of bin ids for each value
	bin_ids = zeros(Int, size(values))

	min, max = extrema(values)
	if min == max
		# If values are all the same, assign them all to bin 1
		number_of_bins = 1
		bin_ids += 1
	elseif mode == "uniform_width"
		bin_ids += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values)), values)
	elseif mode == "uniform_count"
		# TODO: look into why DiscretizeUniformCount doesn't always work, and fix if possible
		try
			bin_ids += encode(LinearDiscretizer(binedges(DiscretizeUniformCount(number_of_bins), values)), values)
		catch
			bin_ids += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values)), values)
			println("Uniform count failed, fell back to uniform width")
		end
	# TODO: include "bayesian_blocks"
	# TODO: Handle "mode doesn't exist" error
	end

	return bin_ids, number_of_bins
end