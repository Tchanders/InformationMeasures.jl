# Tests for partial information decomposition function

# Test partial information decomposition for AND gate
arr_and = [arr1[i] >= 0.5 && arr2[i] >= 0.5 ? 1.0 : 0.0 for i in 1:length(arr1)]
pid = get_partial_information_decomposition(arr1, arr2, arr_and)
@test_approx_eq_eps pid["redundancy"] 0.3 0.05
@test_approx_eq_eps pid["unique_1"] 0 0.05
@test_approx_eq_eps pid["unique_2"] 0 0.05
@test_approx_eq_eps pid["synergy"] 0.5 0.05
println("Partial information decomposition with AND gate passed.")

# Test partial information decomposition for AND gate
pid = get_partial_information_decomposition(arr1, arr2, arr_and, all_orientations = true)
@test_approx_eq_eps pid["z"]["redundancy"] 0.3 0.05
@test_approx_eq_eps pid["z"]["unique_1"] 0 0.05
@test_approx_eq_eps pid["z"]["unique_2"] 0 0.05
@test_approx_eq_eps pid["z"]["synergy"] 0.5 0.05
println("Partial information decomposition with AND gate (all orientations) passed.")