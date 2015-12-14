# Check entropy is log_base(b) for uniform distribution with b bins
@test_approx_eq_eps getentropymillermadow(arr) log2(b) 0.5
@test_approx_eq_eps getentropymillermadow(arr, e) log(b) 0.5

println("Miller Madow entropy passed")