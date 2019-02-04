# Tests for mutual information function
# (NB These tests will fail with the shrinkage estimator, because the 1D,
# 2D and 3D distributions on the right will shrink towards different
# targets by default, whereas the mutual information formula overcomes
# by estimating the probabilities of the 3D distribution, then projecting)

mi = entropy1 + entropy2 - entropy12

# Test mutual information
@test get_mutual_information(arr1, arr2) ≈ mi
println("Mutual information passed")

# Test mutual information formula
@test apply_mutual_information_formula(entropy1, entropy2, entropy12) ≈ mi
println("Mutual information formula passed.")
