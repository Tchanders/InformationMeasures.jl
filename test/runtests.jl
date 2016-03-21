using Base.Test

include("../src/Discretization.jl")
include("../src/Estimators.jl")
include("../src/Formulae.jl")
include("../src/ExportedFunctions.jl")

arr1 = rand(100)
arr2 = rand(100)
arr3 = rand(100)

include("testDiscretization.jl")
include("testEntropy.jl")
include("testConditionalEntropy.jl")
include("testMutualInformation.jl")
include("testConditionalMutualInformation.jl")
include("testInteractionInformation.jl")
include("testTotalCorrelation.jl")
