# Tests concerning cross-entropy

arr21 = rand(100)
arr22 = rand(100)
arr41 = rand(10000)
arr42 = rand(10000)

# The cross-entropy between two copies of an array
@test cross_entropy(arr41, arr41) \approx log2(100) atol=0.01

# The cross-entropy between an array and itself is equal to the regular entropy of that array
@test cross_entropy(arr21, arr21) == get_entropy(arr21)

# The cross-entropy between two small samples from the same distribution will be somewhat larger than the regular entropy
@test cross_entropy(arr21, arr22) \approx log2(10) atol=0.4
@test cross_entropy(arr41, arr42) \approx log2(100) atol=0.02
