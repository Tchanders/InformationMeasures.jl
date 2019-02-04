using InformationMeasures
using Test

arr1 = rand(1000)
arr2 = rand(1000)
arr3 = rand(1000)

# Variables are shared between these tests - don't re-order
include("testDiscretization.jl")
include("testEntropy.jl")
include("testConditionalEntropy.jl")
include("testMutualInformation.jl")
include("testConditionalMutualInformation.jl")
include("testInteractionInformation.jl")
include("testTotalCorrelation.jl")
include("testPartialInformationDecomposition.jl")
include("testCrossEntropy.jl")
