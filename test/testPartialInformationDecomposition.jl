# Tests for partial information decomposition function

# Test partial information decomposition for AND gate
arr_and = [arr1[i] >= 0.5 && arr2[i] >= 0.5 ? 1.0 : 0.0 for i in 1:length(arr1)]
pid = get_partial_information_decomposition(arr1, arr2, arr_and)
@test_approx_eq_eps pid[1] 0.3 0.05
@test_approx_eq_eps pid[2] 0 0.05
@test_approx_eq_eps pid[3] 0 0.05
@test_approx_eq_eps pid[4] 0.5 0.05
println("Partial information decomposition with AND gate passed.")