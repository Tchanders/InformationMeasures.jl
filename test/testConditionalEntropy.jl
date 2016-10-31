# Tests for conditional entropy function

# Test with 1 values array conditioned on itself
@test_approx_eq_eps get_conditional_entropy(arr1, arr1) 0 0.05
println("Conditional entropy with 1 array passed.")

# Test conditional entropy formula
@test_approx_eq_eps get_conditional_entropy(arr1, arr2) get_entropy(arr1, arr2) - get_entropy(arr2) 0.05
println("Conditional entropy formula passed.")
