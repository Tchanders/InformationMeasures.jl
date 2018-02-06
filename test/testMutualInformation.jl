# Tests for mutual information function
# (NB These tests will fail with the shrinkage estimator, because the 1D,
# 2D and 3D distributions on the right will shrink towards different
# targets by default, whereas the mutual information formula overcomes
# by estimating the probabilities of the 3D distribution, then projecting)

# Test mutual information formula
@test get_mutual_information(arr1, arr2) ≈ get_entropy(arr1) + get_entropy(arr2) - get_entropy(arr1, arr2)
println("Mutual information formula passed.")
