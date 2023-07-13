module PhotoAmbiguities

using Parameters
using PartialWaveFunctions
using DataFrames
using Optim
using LinearAlgebra
using StaticArrays

export ε, d, Y
include("specialfunctions.jl")

export Wave, update
export Model
export intensity
include("waves.jl")

export generate
export arg, scale
include("stattools.jl")

export Experiment
export coordinates
export NNL
export fold, unfold
export go2min
export track 
include("experiments.jl")

end # module PhotoAmbiguities
