# Tests for mutual information function
# (NB These tests will fail with the shrinkage estimator, because the 1D,
# 2D and 3D distributions on the right will shrink towards different
# targets by default, whereas the mutual information formula overcomes
# by estimating the probabilities of the 3D distribution, then projecting)

# Test mutual information between two identical arrays is the same as the
# entropy of one of them
@test get_mutual_information(arr1, arr1) == get_entropy(arr1)
println("Mutual information of 2 identical arrays passed.")

# Test mutual information formula
@test get_mutual_information(arr1, arr2) == get_entropy(arr1) + get_entropy(arr2) - get_entropy(arr1, arr2)
println("Mutual information formula passed.")
