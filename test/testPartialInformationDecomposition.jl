# Tests for partial information decomposition function

# Test partial information decomposition for AND gate
arr_and = [arr1[i] >= 0.5 && arr2[i] >= 0.5 ? 1.0 : 0.0 for i in 1:length(arr1)]
pid = get_partial_information_decomposition(arr1, arr2, arr_and)
@test pid["redundancy"] ≈ 0.3 atol = 0.1
@test pid["unique_1"] ≈ 0 atol = 0.1
@test pid["unique_2"] ≈ 0 atol = 0.1
@test pid["synergy"] ≈ 0.5 atol = 0.1
println("Partial information decomposition with AND gate passed.")

pid = get_partial_information_decomposition(arr1, arr2, arr_and, all_orientations = true)
@test pid["z"]["redundancy"] ≈ 0.3 atol = 0.1
@test pid["z"]["unique_1"] ≈ 0 atol = 0.1
@test pid["z"]["unique_2"] ≈ 0 atol = 0.1
@test pid["z"]["synergy"] ≈ 0.5 atol = 0.1
println("Partial information decomposition with AND gate (all orientations) passed.")

# Test redundancy for AND gate
redundancy = get_redundancy(arr1, arr2, arr_and)
@test redundancy ≈ 0.3 atol = 0.1
println("Redundancy with AND gate passed.")

redundancy = get_redundancy(arr1, arr2, arr_and, all_orientations = true)[1]
@test redundancy ≈ 0.3 atol = 0.1
println("Redundancy with AND gate (all orientations) passed.")

# Test unique information for AND gate
mutual_information = get_mutual_information(arr1, arr_and)
@test apply_unique_information_formula(mutual_information, redundancy) ≈ 0 atol = 0.1
println("Unique information with AND gate passed.")

# Test synergy for AND gate
interaction_information = get_interaction_information(arr1, arr2, arr_and)
@test apply_synergy_formula(interaction_information, redundancy) ≈ 0.5 atol = 0.1
println("Synergy with AND gate passed.")
