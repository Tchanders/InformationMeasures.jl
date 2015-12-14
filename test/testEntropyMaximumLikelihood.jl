# Check entropy is log_base(b) for uniform distribution with b bins
@test_approx_eq_eps getentropymaximumlikelihood(arr) log2(b) 0.5
@test_approx_eq_eps getentropymaximumlikelihood(arr, e) log(b) 0.5

println("Maximum likelihood entropy passed")

# # Check mutual information is log_base(b) between two identical
# # variables, uniformly distributed with b bins
# @test_approx_eq_eps mutualinformationmaximumlikelihood(arr, arr) log2(b) 0.5
# @test_approx_eq_eps mutualinformationmaximumlikelihood(arr, arr, e) log(b) 0.5
