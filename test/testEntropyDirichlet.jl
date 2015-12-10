d = 1
n = 100
arr = rand(d, n)

# Check entropy is log_base(b) for uniform distribution with b bins
b = sqrt(n) # getfrequencies function uses sqrt(n) bins
# some choices for a:
# a = 0          :   empirical estimate
# a = 1          :   Laplace
# a = 1/2        :   Jeffreys
# a = 1/m        :   Schurmann-Grassberger  (m: number of bins)
# a = sqrt(n)/m  :   minimax
for a in [0, 1, 0.5, 1/b]
	@test_approx_eq_eps getentropydirichlet(arr, a) log2(b) 0.5
	@test_approx_eq_eps getentropydirichlet(arr, a, e) log(b) 0.5
end

println("Dirichlet passed")
