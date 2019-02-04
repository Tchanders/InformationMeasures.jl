# Tests concerning cross-entropy

arr21 = rand(100)
arr22 = rand(100)
arr41 = rand(10000)
arr42 = rand(10000)

# The cross-entropy between two copies of an array
@test get_cross_entropy(arr41, arr41) ≈ log2(100) atol = 0.01
println("Cross-entropy between array and itself passed.")

# The cross-entropy between an array and itself is equal to the regular entropy of that array
@test get_cross_entropy(arr21, arr21) == get_entropy(arr21)
println("Equality of cross-entropy between array and itself and entropy passed.")

# The cross-entropy between two small samples from the same distribution will be somewhat larger than the regular entropy
@test get_cross_entropy(arr21, arr22) ≈ log2(10) atol = 0.4
println("Cross-entropy between two small sample arrays passed.")
@test get_cross_entropy(arr41, arr42) ≈ log2(100) atol = 0.02
println("Cross-entropy between two large sample arrays passed.")
