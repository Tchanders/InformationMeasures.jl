d = 1
n = 100
arr = rand(d, n)

# Check DomainError is thrown when !(k < n)
@test_throws DomainError entropyknn(arr, 100)

# Check no DomainError is thrown when k < n
@test typeof(entropyknn(arr, 99)) == Float64

# Check entropy is log_base(n) for uniform distribution width n
for m in [10 100 1000]
	arr = rand(d, n) * m
	@test_approx_eq_eps entropyknn(arr) log2(m) 0.5
	@test_approx_eq_eps entropyknn(arr, 3, 10) log10(m) 0.5
	@test_approx_eq_eps entropyknn(arr, 3, e) log(m) 0.5
end
