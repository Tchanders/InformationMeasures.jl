# InformationMeasures

*Release version:*

[![InformationMeasures](http://pkg.julialang.org/badges/InformationMeasures_0.4.svg)](http://pkg.julialang.org/?pkg=InformationMeasures)
[![InformationMeasures](http://pkg.julialang.org/badges/InformationMeasures_0.5.svg)](http://pkg.julialang.org/?pkg=InformationMeasures)
[![InformationMeasures](http://pkg.julialang.org/badges/InformationMeasures_0.6.svg)](http://pkg.julialang.org/?pkg=InformationMeasures)

*Development version:*

[![Build Status](https://travis-ci.org/Tchanders/InformationMeasures.jl.svg?branch=master)](https://travis-ci.org/Tchanders/InformationMeasures.jl)
[![codecov.io](http://codecov.io/github/Tchanders/InformationMeasures.jl/coverage.svg?branch=master)](http://codecov.io/github/Tchanders/InformationMeasures.jl?branch=master)

## Installation

`Pkg.add("InformationMeasures")`

## Performance

In cases where optimal performance is needed, the latest version of InformationMeasures is recommended, with Julia 1.0 or later. See also [Advanced usage](#advanced-usage).

## Basic usage

Currently information measures on three or fewer variables are supported. The basic use case is to pass data arrays for each variable into each function. These will be discretized.

It is also possible to pass in frequencies (if the data has already been discretized), or probabilities (if the probabilities are already known or have already been estimated) - see below.

```
using InformationMeasures

data_1 = rand(100)
data_2 = rand(100)
data_3 = rand(100)

# Entropy
ent_1 = get_entropy(data_1)
ent_12 = get_entropy(data_1, data_2)
ent_123 = get_entropy(data_1, data_2, data_3)

# Conditional entropy
ce_1_on_2 = get_conditional_entropy(data_1, data_2)

# Mutual information
mi_12 = get_mutual_information(data_1, data_2)

# Conditional mutual information
cmi_12_on_3 = get_conditional_mutual_information(data_1, data_2, data_3)

# Interaction information
ii_123 = get_interaction_information(data_1, data_2, data_3)

# Total correlation
tc_123 = get_total_correlation(data_1, data_2, data_3)

# Partial information decomposition
pid_123 = get_partial_information_decomposition(data_1, data_2, data_3)
```

## Config options

The following keyword arguments can be passed in to each function:

**estimator** (String) Estimator for estimating the probability distribution
* `"maximum_likelihood"` (default)
* `"miller_madow"`
* `"dirichlet"`
* `"shrinkage"`

**base** (Number) Base of the logarithm, i.e. the units for entropy
* `2` (default)

**mode** (String) Method for discretizing
* `"uniform_width"` (default)
* `"uniform_count"`
* `"bayesian_blocks"`

**number_of_bins** (Integer)
* `0` (default)

**get_number_of_bins** (Function) Customized function for calculating the number of bins (will only be used if `number_of_bins` is `0`)
* `get_root_n` (default)

#### Estimator-specific config options

**lambda** (Void or Number) Shrinkage intensity (if left as `nothing`, will be calculated automatically)
* `nothing` (default)

**prior** (Number) Dirichlet prior (if left as `0`, Dirichlet estimator is equivalent to maximum likelihood)
* `0` (default)

#### Values, frequencies, or probabilities

The information measures can be calculated from raw data values, frequencies (if the data has already been discretized), or probabilities (if the probabilities are already known or have already been estimated).

To calculate entropy from frequencies, call `get_entropy` with the keyword argument `discretized = true`

For all other information measures, simply pass in a single array of frequencies or probabilities (2D for conditional entropy and mutual information or 3D for conditional mutual information, mutual information and total correlation). If they are probabilities, include the keyword argument `probabilities = true`, otherwise they will be treated as frequencies.

## Discretization

Although discretization is taken care of when the information measures are calculated, it is possible to discretize raw values directly, for example to investigate how different discretization algorithms and bin numbers affect the discretization.

```
data = rand(100)
disc_val = discretize_values(data)
```

NB `discretize_values` returns the frequencies for each bin in order, rather than the discretized values. An example of how to get the discretized values is discussed below.

## Advanced usage

Functions such as `get_entropy` and `get_mutual_information` are designed to be flexible and easy to use with different types of input and config options. In some cases it may be quicker to bypass these functions.

### Example

When calculating the mutual information between every pair of data vectors from a large dataset, simply calling `get_mutual_information` on each pair of vectors will result in each vector being discretized multiple times.

Currently, discretization for multiple variables works by discretizing the marginals independently, then reconstructing the higher dimensional frequencies from these discretized marginals. Therefore discretizing each variable once in advance will not affect the results, but will be much quicker. Joint frequencies cannot be reconstructed from the bin frequencies; instead the discretized values should be stored. `get_bin_ids!` should therefore be used, instead of `discretize_values`:

```
data_1 = rand(100)
data_2 = rand(100)
data_3 = rand(100)

number_of_bins = 10
mi_base = 2

bin_ids_1 = zeros(Int, length(data_1))
get_bin_ids!(data_1, "uniform_width", number_of_bins, bin_ids_1)

bin_ids_2 = zeros(Int, length(data_2))
get_bin_ids!(data_2, "uniform_width", number_of_bins, bin_ids_2)

bin_ids_3 = zeros(Int, length(data_3))
get_bin_ids!(data_3, "uniform_width", number_of_bins, bin_ids_3)

f_12 = get_frequencies_from_bin_ids(bin_ids_1, bin_ids_2, number_of_bins, number_of_bins)
p_12 = get_probabilities("maximum_likelihood", f_12)
mi_12 = apply_mutual_information_formula(p_12, sum(p_12, dims = 2), sum(p_12, dims = 1), mi_base)

f_13 = get_frequencies_from_bin_ids(bin_ids_1, bin_ids_3, number_of_bins, number_of_bins)
p_13 = get_probabilities("maximum_likelihood", f_13)
mi_13 = apply_mutual_information_formula(p_13, sum(p_13, dims = 2), sum(p_13, dims = 1), mi_base)

# And so on...
```

Note that the probability distribution is estimated from the joint frequencies rather than the marginals, meaning that, for most estimators, `sum(p_12, dims = 2)` may be different from `sum(p_13, dims = 2)`, even though both represent the estimated probability distribution for `data_1`. (This is not the case for the "maximum likelihood" estimator, which just divides the bin frequencies by the total frequency. For this estimator, the marginal probabilities could be stored in advance to avoid calculating them as they are passed into `apply_entropy_formula`. The best performance in that case may depend on the cost of storage vs calculations.)

Here are two full examples of the "quick" vs the "easy" way to estimate the mutual information between all pairs of a set of variables.

```
data = rand(100, 100)

function mi_quick(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2)

	nvals, nvars = size(data)

	bin_ids = zeros(Int, (nvals, nvars))
	nbins = Int(round(sqrt(nvals)))
	mis = zeros(binomial(nvars, 2))

	for i in 1 : nvars
		get_bin_ids!(data[1:nvals, i:i], discretizer, nbins, view(bin_ids, 1:nvals, i:i))
	end

	index = 1
	for i in 1 : nvars, j in i+1 : nvars
		f = get_frequencies_from_bin_ids(bin_ids[1:end, i:i], bin_ids[1:end, j:j], nbins, nbins)
		p = get_probabilities(estimator, f)
		mis[index] = apply_mutual_information_formula(p, sum(p, dims = 1), sum(p, dims = 2), mi_base)
		index += 1
	end

	return mis

end

function mi_easy(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2)
	nvals, nvars = size(data)
	mis = zeros(binomial(nvars, 2))

	index = 1
	for i in 1 : nvars, j in i+1 : nvars
		mis[index] = get_mutual_information(data[1:nvals, i:i], data[1:nvals, j:j]; mode = discretizer, estimator = estimator, base = mi_base)
		index += 1
	end

	return mis
end
```

## Contributing

Contributions and bug reports are welcome!

## References

[1] Chan, Stumpf and Babtie (2017) [Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures](http://www.cell.com/cell-systems/fulltext/S2405-4712(17)30386-1) Cell Systems
