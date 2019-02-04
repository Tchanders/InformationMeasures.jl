# Tests for interaction information function

ii = cmi - mi

# Test interaction information
@test get_interaction_information(arr1, arr2, arr3) ≈ ii
println("Interaction information passed.")

# Test interaction information formula
@test apply_interaction_information_formula(cmi, mi) ≈ ii
println("Interaction information formula passed.")
