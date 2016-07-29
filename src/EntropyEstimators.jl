module EntropyEstimators

# The Bayesian blocks function should be incorporated into Discretizers.jl
# For now it is inluded here separately
include("../src/BayesianBlocks.jl")

include("../src/Discretization.jl")
include("../src/Estimators.jl")
include("../src/Formulae.jl")
include("../src/ExportedFunctions.jl")

end # module
