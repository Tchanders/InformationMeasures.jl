# Tests for conditional entropy function

entropy1 = get_entropy(arr1)
entropy2 = get_entropy(arr2)
entropy12 = get_entropy(arr1, arr2)

# Test with 1 values array conditioned on itself
@test get_conditional_entropy(arr1, arr1) ≈ 0 atol = 0.05
println("Conditional entropy with 1 array passed.")

conditional_entropy = entropy12 - entropy2

# Test conditional entropy
@test get_conditional_entropy(arr1, arr2) ≈ conditional_entropy atol = 0.05
println("Conditional entropy passed.")

# Test conditional entropy formula
@test apply_conditional_entropy_formula(entropy12, entropy2) ≈ conditional_entropy atol = 0.05
println("Conditional entropy formula passed.")
