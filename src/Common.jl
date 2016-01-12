# WIP The functions common to different estimators

"""
Coverts observed values into bin IDs, for discretization into root(n)
discrete bins, where n is the number of observed values.

Parameters:

values - Array{Float64,2} - The observed values.
"""
function getbinids(values::Array{Float64,2})

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


"""
Plugs probabilities into the entropy formula
"""
function applyentropyformula(probabilities, base)
	return -sum([p == 0 ? 0 : p .* log(p) for p in probabilities]) / log(base)
end
