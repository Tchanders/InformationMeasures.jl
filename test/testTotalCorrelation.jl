# Tests for the total correlation function

total_correlation = entropy1 + entropy2 + entropy3 - entropy123

# Test total correlation
@test get_total_correlation(arr1, arr2, arr3) ≈ total_correlation
println("Total correlation passed.")

# Test total correlation formula
@test apply_total_correlation_formula(entropy1, entropy2, entropy3, entropy123) ≈ total_correlation
println("Total correlation formula passed.")
