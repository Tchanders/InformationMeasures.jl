# Test discretization functions

arr = rand(100)

# Test with 1 values array
@test size(discretize_values(arr)) == (10, 1)
println("Discretize values with 1 array passed.")

# Test with 2 values arrays
@test size(discretize_values(arr, arr)) == (10, 10)
println("Discretize values with 2 arrays passed.")

# Test with 3 values arrays
@test size(discretize_values(arr, arr, arr)) == (10, 10, 10)
println("Discretize values with 3 arrays passed.")

# Test with 4 values arrays
@test size(discretize_values(arr, arr, arr, arr)) == (10, 10, 10, 10)
println("Discretize values with >3 arrays passed.")

# Test with uniform_width
@test size(discretize_values(arr, mode = "uniform_width")) == (10, 1)
println("Discretize values with uniform width passed.")

# Test with uniform_count
@test size(discretize_values(arr, mode = "uniform_count")) == (10, 1)
println("Discretize values with uniform count passed.")

# Test with bayesian_blocks
# @test size(discretize_values(arr, mode = "bayesian_blocks")) == (10, 1)
# println("Discretize values with Bayesian blocks passed.")

# Test with number_of_bins
@test size(discretize_values(arr, number_of_bins = 5)) == (5, 1)
println("Discretize values with number of bins passed.")

# Test with get_number_of_bins
@test size(discretize_values(arr, get_number_of_bins = function(values) return 2 end)) == (2, 1)
println("Discretize values with get number of bins passed.")
