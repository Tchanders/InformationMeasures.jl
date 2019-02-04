# Tests for the conditional mutual information function

entropy3 = get_entropy(arr3)
entropy13 = get_entropy(arr1, arr3)
entropy23 = get_entropy(arr2, arr3)
entropy123 = get_entropy(arr1, arr2, arr3)
cmi = entropy13 + entropy23 - entropy123 - entropy3

# Test conditional mutual information
@test get_conditional_mutual_information(arr1, arr2, arr3) ≈ cmi
println("Conditional mutual information passed.")

# Test conditional mutual information formula
@test apply_conditional_mutual_information_formula(entropy13, entropy23, entropy123, entropy3) ≈ cmi
println("Conditional mutual information formula passed.")
