# Estimates the entropy of a continuous data set using the
# algorithm outlined in: Kraskov A, Stögbauer H, Grassberger P.
# Estimating mutual information. Physical Review E. 2004;
# 69(6):066138, previously implemented in Python:
# http://www.isi.edu/~gregv/npeet.html.

export getentropyknn

using NearestNeighbors

"""
Calculates an estimate for the entropy of a set of observed
values using a k-nearest neighbour algorithm. First noise is
added to the values, then a k-d tree is made, then the tree
is queried for the distance to the kth nearest neighbour of
each value, then these distances are substituted into the
formula.

Parameters:

values - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions of each value and n is the number of
values. (The KDTree module requires an array of this shape.)

k - Int - The degree of nearest neighbours to find. Must be
smaller than n.

base - Int - The base of the logarithm, i.e. the units.

intens - Float - A multiplier for adding noise.
"""
function getentropyknn(values::Array{Float64,2}, k=3, base=2, noise=1e-10)
	dimensions, n = size(values)
	assert(k < n)

	# Add noise to each value
	values = values + rand(dimensions, n) * noise

	# Make the KD tree
	tree = KDTree(values)

	# Make array of distances to kth nearest neighbour for each
	# value
	distances = Array(Float64, 0)
	for i in 1:n
		# Query the tree for the k nearest neighbours using knn
		# (requires the values array to be reshaped)
		point = reshape(values[1:dimensions, i:i], dimensions)
		push!(distances, knn(tree, point, k + 1, true)[2][k + 1])
	end

	# Substitute distances into the equation
	constant = digamma(n) - digamma(k) + dimensions * log(2)
	return (constant + dimensions * mean(log(distances))) / log(base)
end
