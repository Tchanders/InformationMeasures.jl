# Tests for interaction information function

# Test interaction information formula
@test_approx_eq get_interaction_information(arr1, arr2, arr3) get_conditional_mutual_information(arr1, arr2, arr3) - get_mutual_information(arr1, arr2)
println("Interaction information formula passed.")