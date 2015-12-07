# The functions common to different estimators

# Needs improving
function discretizecounts(counts::Array{Float64,2})
	x = reshape(counts, length(counts))
	numBins = round(sqrt(length(x)))
	sort!(x)
	min, max = extrema(x)
	diff = (max - min) / numBins
	r = Range(min : diff : max)
	c = zeros(length(r) - 1)
	for i in 1 : length(c)
		while length(x) > 0 && x[1] <= r[i + 1]
			c[i] += 1
			shift!(x)
		end
	end
	# Range sometimes doesn't go up to max. Hack to add an extra bin if so, and
	# add any leftovers to the last bin:
	# If there are leftovers...
	if length(x) > 0
		# ...if there are the correct number of bins...
		if length(c) === numBins
			# ...then add the leftovers to the last bin...
			c[end] += length(x)
		# ...but if we are one bin short...
		elseif length(c) === numBins - 1
			# ...add a new bin with the number left over.
			push!(c, length(x))
		end
	end
	c
end

function entropyformula(frequencies, base)
	return -sum(frequencies .* log(frequencies)) / log(base)
end
