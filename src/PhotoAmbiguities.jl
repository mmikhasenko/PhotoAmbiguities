module PhotoAmbiguities

using Parameters
using PartialWaveFunctions
using DataFrames
using Optim
using LinearAlgebra
using StaticArrays
using Interpolations
using ForwardDiff

export ε, d, Y
include("specialfunctions.jl")

export Wave, update
export Model
export intensity
export standardize
include("waves.jl")

export generate
export arg, scale
export distance
export cluster
export add_conjugate
export drop_conjugate
include("stattools.jl")

export Experiment
export coordinates
export NNL
export fold, unfold
export go2min, go2min_precompute
export nminima
export pullbilinear
export BiLinearCompute
include("experiments.jl")

end # module PhotoAmbiguities
