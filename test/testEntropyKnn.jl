# Check DomainError is thrown when !(k < n)
@test_throws AssertionError getentropyknn(arr, 100)

# Check no DomainError is thrown when k < n
@test typeof(getentropyknn(arr, 99)) == Float64

# Check entropy is log_base(n) for uniform distribution width n
for m in [10 100 1000]
	arr = rand(d, n) * m
	@test_approx_eq_eps getentropyknn(arr) log2(m) 0.5
	@test_approx_eq_eps getentropyknn(arr, 3, e) log(m) 0.5
end

println("K-nearest neighbor entropy passed")
