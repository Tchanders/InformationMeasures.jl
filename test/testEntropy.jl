# Tests for entropy function

arr = rand(1000)
sqrt_1000 = sqrt(1000)

# Test with 1 values array
@test_approx_eq_eps get_entropy(arr) log2(sqrt_1000) 0.05
println("Entropy with 1 array passed.")

# Test with 2 values arrays
@test_approx_eq_eps get_entropy(arr, arr) log2(sqrt_1000) 0.05
println("Entropy with 2 arrays passed.")

# Test with 3 values arrays
@test_approx_eq_eps get_entropy(arr, arr, arr) log2(sqrt_1000) 0.05
println("Entropy with 3 arrays passed.")

# Test with uniform_width
@test_approx_eq_eps get_entropy(arr, mode = "uniform_width") log2(sqrt_1000) 0.05
println("Entropy with uniform width passed.")

# Test with uniform_count
@test_approx_eq_eps get_entropy(arr, mode = "uniform_count") log2(sqrt_1000) 0.05
println("Entropy with uniform count passed.")

# Test with bayesian_blocks
# @test_approx_eq_eps get_entropy(arr, mode = "bayesian_blocks") log2(10) 0.05
# println("Entropy with Bayesian blocks passed.")

# Test with number_of_bins
@test_approx_eq_eps get_entropy(arr, number_of_bins = 5) log2(5) 0.05
println("Entropy with number of bins passed.")

# Test with get_number_of_bins
@test_approx_eq_eps get_entropy(arr, get_number_of_bins = function(values) return 2 end) log2(2) 0.05
println("Entropy with get number of bins passed.")

# Test with maximum_likelihood
@test_approx_eq_eps get_entropy(arr, estimator = "maximum_likelihood") log2(sqrt_1000) 0.05
println("Entropy with maximum likelihood passed.")

# Test with shrinkage
@test_approx_eq_eps get_entropy(arr, estimator = "shrinkage") log2(sqrt_1000) 0.05
println("Entropy with shrinkage passed.")

# Test with shrinkage and lambda
@test_approx_eq_eps get_entropy(arr, estimator = "shrinkage", lambda = 0) log2(sqrt_1000) 0.05
println("Entropy with shrinkage and lambda passed.")

# Test with dirichlet
@test_approx_eq_eps get_entropy(arr, estimator = "dirichlet") log2(sqrt_1000) 0.05
println("Entropy with Dirichlet passed.")

# Test with dirichlet and prior
@test_approx_eq_eps get_entropy(arr, estimator = "dirichlet", prior = 1) log2(sqrt_1000) 0.05
println("Entropy with Dirichlet and prior passed.")

# Test with miller_madow
@test_approx_eq_eps get_entropy(arr, estimator = "miller_madow") log2(sqrt_1000) 0.05
println("Entropy with Miller-Madow passed.")

# Test with base e
@test_approx_eq_eps get_entropy(arr, base = e) log(sqrt_1000) 0.05
println("Entropy with change of base passed.")

# Test with discretized
@test_approx_eq_eps get_entropy(discretize_values(arr), discretized = true) log2(sqrt_1000) 0.05
println("Entropy with discretized input passed.")
