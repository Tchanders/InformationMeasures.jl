# Tests for the conditional mutual information function

# Test conditional mutual information formula
@test get_conditional_mutual_information(arr1, arr2, arr3) ≈ get_entropy(arr1, arr3) + get_entropy(arr2, arr3) - get_entropy(arr1, arr2, arr3) - get_entropy(arr3)
println("Conditional mutual information formula passed.")
