module PhotoAmbiguities

using Parameters
using PartialWaveFunctions
using DataFrames
using Optim
using LinearAlgebra
using StaticArrays
using Interpolations

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
export cluster, drop_conjugate
include("stattools.jl")

export Experiment
export coordinates
export NNL
export fold, unfold
export go2min
export nminima
include("experiments.jl")

end # module PhotoAmbiguities
