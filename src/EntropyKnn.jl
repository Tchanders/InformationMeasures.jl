# Estimates the entropy of a continuous data set using the
# algorithm outlined in: Kraskov A, Stögbauer H, Grassberger P.
# Estimating mutual information. Physical Review E. 2004;
# 69(6):066138, previously implemented in Python:
# http://www.isi.edu/~gregv/npeet.html.

export getentropyknn

using NearestNeighbors

"""
Estimates the entropy of a continuous data set using a
K-nearest neighbour algorithm.

Parameters:

data - dxn Array{Float64,2} - A data array where n is the number
of data points and d is the dimensionality of the data. KDTree
requires an array of this shape.

k - Int - The degree of nearest neighbours to find. Must be
smaller n.

base - Int - The base of the logarithm, i.e. the units.

intens - Float - A multiplier for adding noise.
"""
function getentropyknn(data::Array{Float64,2}, k=3, base=2, noise=1e-10)
	dimensions, n = size(data)
	assert(k < n)

	# Add noise to each data point
	data = data + rand(dimensions, n) * noise

	# Make the KD tree
	tree = KDTree(data)

	# Make array of distances to kth nearest neighbour for each
	# data point
	distances = Array(Float64, 0)
	for i in 1:n
		# Query the tree for the k nearest neighbours using knn
		# (requires the data to be reshaped)
		point = reshape(data[1:dimensions, i:i], dimensions)
		push!(distances, knn(tree, point, k + 1, true)[2][k + 1])
	end

	# Substitute distances into the equation
	constant = digamma(n) - digamma(k) + dimensions * log(2)
	return (constant + dimensions * mean(log(distances))) / log(base)
end
