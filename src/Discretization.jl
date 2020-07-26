# Functions for discretizing continuous values using various algorithms

using Discretizers

export get_bin_ids!, get_frequencies, get_frequencies_from_bin_ids

# Parameters:
# 	- values_x, arrays of floats (multiple arrays supported)
function get_root_n(values_x...)
	return round(Int, sqrt(size(values_x[1])[1]))
end

# Parameters:
# 	- values_x, arrays of floats (multiple arrays supported)
# 	- mode, string
#	- number_of_bins, number
function get_frequencies(mode, number_of_bins, values_x)
	bin_ids = zeros(Int, size(values_x))
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
	bin_ids_x = zeros(Int, size(values_x))
	bin_ids_y = zeros(Int, size(values_y))
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
	bin_ids_x = zeros(Int, size(values_x))
	bin_ids_y = zeros(Int, size(values_y))
	bin_ids_z = zeros(Int, size(values_z))
	number_of_bins_x = get_bin_ids!(values_x, mode, number_of_bins, bin_ids_x)
	number_of_bins_y = get_bin_ids!(values_y, mode, number_of_bins, bin_ids_y)
	number_of_bins_z = get_bin_ids!(values_z, mode, number_of_bins, bin_ids_z)
	frequencies = zeros(Int, (number_of_bins_x, number_of_bins_y, number_of_bins_z))
	for i in 1:n
		frequencies[bin_ids_x[i], bin_ids_y[i], bin_ids_z[i]] += 1
	end
	return frequencies
end
function get_frequencies(mode, number_of_bins, values_x...)
	# n is the number of data points
	n = size(values_x[1])[1]
	# d is the number of dimensions
	d = length(values_x)
	bin_ids = [zeros(Int, 1) for i in 1:d]
	all_number_of_bins = zeros(Int, d)
	for i in 1:d
		current_values = values_x[i]
		bin_ids[i] = zeros(Int, size(current_values))
		all_number_of_bins[i] = get_bin_ids!(current_values, mode, number_of_bins, bin_ids[i])
	end
	frequencies = zeros(Int, tuple(all_number_of_bins...))
	for i in 1:n
		frequencies[ntuple(x -> bin_ids[x][i], d)...] += 1
	end
	return frequencies
end

# Parameters:
# 	- bin_ids_x, array of ints
# 	- number_of_bins_x, number
function get_frequencies_from_bin_ids(bin_ids_x, number_of_bins_x)
	n = size(bin_ids_x)[1]
	frequencies = zeros(Int, number_of_bins_x)
	for i in 1:n
		frequencies[bin_ids_x[i]] += 1
	end
	return frequencies
end
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
# 	- values_x, array of floats
# 	- mode, string
#	- number_of_bins, integer
#	- bin_ids, 1-dimensional array of bin ids for each value
function get_bin_ids!(values_x, mode, number_of_bins, bin_ids)
	min, max = extrema(values_x)
	if min == max
		# If values are all the same, assign them all to bin 1
		number_of_bins = 1
		bin_ids[1:end] .= 1
	elseif mode == "uniform_width"
		bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values_x)), values_x)
	elseif mode == "binarize"
		number_of_bins = 2
		bin_ids[1:end] = map(v -> v == 0 ? 1 : 2, values_x)
	elseif mode == "uniform_count"
		try
			bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformCount(number_of_bins), reshape(values_x, length(values_x)))), values_x)
		catch
			bin_ids = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values_x)), values_x)
			println("Uniform count failed, fell back to uniform width")
		end
	elseif mode == "bayesian_blocks"
		try
			edges = binedges(DiscretizeBayesianBlocks(), values_x)
			bin_ids[1:end] = encode(LinearDiscretizer(edges), values_x)
			number_of_bins = length(edges) - 1
		catch
			bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values_x)), values_x)
			println("Bayesian blocks failed, fell back to uniform width")
		end
	else
		println("Mode doesn't exist, fell back to uniform width")
		bin_ids[1:end] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(number_of_bins), values_x)), values_x)
	end

	return number_of_bins
end
