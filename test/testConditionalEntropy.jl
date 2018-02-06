# Tests for conditional entropy function

# Test with 1 values array conditioned on itself
@test get_conditional_entropy(arr1, arr1) ≈ 0 atol = 0.05
println("Conditional entropy with 1 array passed.")

# Test conditional entropy formula
@test get_conditional_entropy(arr1, arr2) ≈ get_entropy(arr1, arr2) - get_entropy(arr2) atol = 0.05
println("Conditional entropy formula passed.")
