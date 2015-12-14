# The functions common to different estimators

# Needs improving
"""
Sorts counts in to frequencies for root(n) bins
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

"""
Plugs probabilities into the entropy formula
"""
function applyentropyformula(probabilities, base)
	return -sum(probabilities .* log(probabilities)) / log(base)
end

##############################

function applymutualinformationformula(freqs2d, base=2)
# 	freqs2d = freqs2d / sum(freqs2d) # check this
	freqsX = sum(freqs2d, 2)
	freqsY = sum(freqs2d, 1)
	freqsNull = (freqsX * freqsY) # Why not transpose?
	return kldivergenceformula(freqs2d, freqsNull, base)
end

function kldivergenceformula(freqs1, freqs2, base=2)

	function zeroOrLess(a)
		return a <= 0
	end

	freqs1 = freqs1 / sum(freqs1)
	freqs2 = freqs2 / sum(freqs2)

	if findfirst(zeroOrLess, freqs2) !== 0
		println("Warning: vanishing values in argument freqs2") # Read up on this
	end

	likelihoodratio = [freqs1[i] > 0 ? log(freqs1[i] / freqs2[i]) : 0 for i = 1:length(freqs1)]
	return sum(freqs1[:] .* likelihoodratio) / log(base)

end
