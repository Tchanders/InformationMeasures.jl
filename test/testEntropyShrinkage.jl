d = 1
n = 10000
arr = rand(d, n)

# Check entropy is log_base(b) for uniform distribution with b bins
b = sqrt(n) # getfrequencies function uses sqrt(n) bins
@test_approx_eq_eps getentropyshrinkage(arr) log2(b) 0.5
@test_approx_eq_eps getentropyshrinkage(arr, e) log(b) 0.5

println("Shrinkage passed")