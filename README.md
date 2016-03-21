# EntropyEstimators

## Installation

Pkg.clone("git://github.com/Tchanders/EntropyEstimators.jl.git")

## Basic usage

Currently information measures on three or fewer variables are supported. The basic use case is to pass data arrays for each variable into each function. These will be discretized.

It is also possible to pass in frequencies (if the data has already been discretized), or probabilities (if the probabilities are already known or have already been estimated) - see below.

```
using EntropyEstimators

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
```

## Config options

The following keyword arguments can be passed in to each function:

**estimator** (String) Estimator for estimating the probability distribution
* "maximum_likelihood" (default)
* "miller_madow"
* "dirichlet"
* "shrinkage"

**base** (Number) Base of the logarithm, i.e. the units for entropy
* 2 (default)

**mode** (String) Method for discretizing
* "uniform_width" (default)
* "uniform_count"

**number_of_bins** (Integer)
* 0 (default)

**get_number_of_bins** (Function) Customized function for calculating the number of bins (will only be used if number_of_bins is 0)
* get_root_n (default)

#### Estimator-specific config options

**lambda** (Void or Number) Shrinkage intensity (if left as `nothing`, will be calculated automatically)
* nothing (default)

**prior** (Number) Dirichlet prior (if left as `0`, Dirichlet estimator is equivalent to maximum likelihood)
* 0 (default)

#### Values, frequencies, or probabilities

The information measures can be calculated from raw data values, frequencies (if the data has already been discretized), or probabilities (if the probabilities are already known or have already been estimated).

To calculate entropy from frequencies, call `get_entropy` with the keyword argument `discretized = true`

For all other information measures, simply pass in a single array of frequencies or probabilities (2D for conditional entropy and mutual information or 3D for conditional mutual information, mutual information and total correlation). If they are probabilities, include the keyword argument `probabilities = true`, otherwise they will be treated as frequencies.

## Discretization

Although discretization is taken care of when the information measures are calculated, it is possible to discretize raw values directly, for example to investigate how different modes and bin numbers affect the discretization.

```
data = rand(100)
disc_val = discretize_values(data)
```
