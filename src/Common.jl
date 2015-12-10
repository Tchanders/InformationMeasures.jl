# The functions common to different estimators

# Needs improving
"""
Gets frequencies from counts, based on root(n) bins
"""
function getfrequencies(counts::Array{Float64,2})
	counts = reshape(counts, length(counts))
	numberofbins = round(sqrt(length(counts)))
	sort!(counts)
	min, max = extrema(counts)
	binwidth = (max - min) / numberofbins
	range = Range(min : binwidth : max)
	frequencies = zeros(length(range) - 1)
	for i in 1 : length(frequencies)
		while length(counts) > 0 && counts[1] <= range[i + 1]
			frequencies[i] += 1
			shift!(counts)
		end
	end
	# Range sometimes doesn't go up to max. Hack to add an extra bin if so, and
	# add any leftovers to the last bin:
	# If there are leftovers...
	if length(counts) > 0
		# ...if there are the correct number of bins...
		if length(frequencies) === numberofbins
			# ...then add the leftovers to the last bin...
			frequencies[end] += length(counts)
		# ...but if we are one bin short...
		elseif length(frequencies) === numberofbins - 1
			# ...add a new bin with the number left over.
			push!(frequencies, length(counts))
		end
	end
	return frequencies
end

function applyentropyformula(frequencies, base)
	return -sum(frequencies .* log(frequencies)) / log(base)
end
