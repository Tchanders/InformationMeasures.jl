# WIP The functions common to different estimators

"""
Coverts observed values into bin IDs, for discretization into root(n)
discrete bins, where n is the number of observed values.

Parameters:

values - Array{Float64,2} - The observed values.
"""
function getbinids(values::Array{Float64,2}) # Why does this sometimes return NaN?

	# This is a separate function because might want to offer different methods in the future
	function getnumberofbins(values)
		# size(values)[2] is n
		return round(Int, sqrt(size(values)[2]))
	end

	# This is a separate function because might want to offer different methods in the future
	function getvalueadjustmentparameters(values, numberofbins)
		min, max = extrema(values)
		binwidth = (max - min) / (numberofbins - 1)
		return min, binwidth
	end

	function adjustvalues(values, numberofbins)
		min, binwidth = getvalueadjustmentparameters(values, numberofbins)
		return floor(Int, (values - min) * (1 / binwidth)) + 1
	end

	numberofbins = getnumberofbins(values)
	numberofdimensions, n = size(values)
	binids = zeros(Int, size(values))
	for i in 1:numberofdimensions
		binids[i:i, 1:end] += adjustvalues(values, numberofbins)
	end
	return binids, numberofbins

end

"""
Converts a set of observed values into frequencies for root(n)
discrete bins, where n is the number of observed values.

Parameters:

values - Array{Float64,2} - The observed values.
"""
function getfrequencies(values::Array{Float64,2})
	binids, numberofbins = getbinids(values)
	frequencies = zeros(Int, (1, numberofbins))
	for binid in binids
		frequencies[binid] += 1
	end
	return frequencies
end

	# Do once for each dimension
	#for i in 1:size(values)[1]

	# NB couns HAS BEEN RENAMED TO values 
	# counts = reshape(counts, length(counts))
	# numberofbins = round(sqrt(length(counts)))
	# sort!(counts)
	# min, max = extrema(counts)
	# binwidth = (max - min) / numberofbins
	# range = Range(min : binwidth : max)
	# frequencies = zeros(length(range) - 1)
	# for i in 1 : length(frequencies)
	# 	while length(counts) > 0 && counts[1] <= range[i + 1]
	# 		frequencies[i] += 1
	# 		shift!(counts)
	# 	end
	# end
	# # Range sometimes doesn't go up to max. Hack to add an extra bin if so, and
	# # add any leftovers to the last bin:
	# # If there are leftovers...
	# if length(counts) > 0
	# 	# ...if there are the correct number of bins...
	# 	if length(frequencies) === numberofbins
	# 		# ...then add the leftovers to the last bin...
	# 		frequencies[end] += length(counts)
	# 	# ...but if we are one bin short...
	# 	elseif length(frequencies) === numberofbins - 1
	# 		# ...add a new bin with the number left over.
	# 		push!(frequencies, length(counts))
	# 	end
	# end
	# return frequencies
#end

# """
# """
# function makejointfrequencies(frequencies1::Array{Float64,1}, frequencies2::Array{Float64,1})
# 	n = length(frequencies1)
# 	frequencymatrix = zeros(n, n)
# 	for i in 1:n
# 		frequencymatrix[frequencies2[i], frequencies1[i]] += 1
# 	end
# 	return frequencymatrix
# end

"""
Plugs probabilities into the entropy formula
"""
function applyentropyformula(probabilities, base)
	return -sum(probabilities .* log(probabilities)) / log(base)
end

# function applymutualinformationformula(freqs2d, base=2)
# # 	freqs2d = freqs2d / sum(freqs2d) # check this
# 	freqsX = sum(freqs2d, 2)
# 	freqsY = sum(freqs2d, 1)
# 	freqsNull = (freqsX * freqsY) # Why not transpose?
# 	return kldivergenceformula(freqs2d, freqsNull, base)
# end

# function kldivergenceformula(freqs1, freqs2, base=2)

# 	function zeroOrLess(a)
# 		return a <= 0
# 	end

# 	freqs1 = freqs1 / sum(freqs1)
# 	freqs2 = freqs2 / sum(freqs2)

# 	if findfirst(zeroOrLess, freqs2) !== 0
# 		println("Warning: vanishing values in argument freqs2") # Read up on this
# 	end

# 	likelihoodratio = [freqs1[i] > 0 ? log(freqs1[i] / freqs2[i]) : 0 for i = 1:length(freqs1)]
# 	return sum(freqs1[:] .* likelihoodratio) / log(base)

# end
