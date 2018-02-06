# Tests for the total correlation function

# Test total correlation formula
@test get_total_correlation(arr1, arr2, arr3) ≈ get_entropy(arr1) + get_entropy(arr2) + get_entropy(arr3) - get_entropy(arr1, arr2, arr3)
println("Total correlation formula passed.")
