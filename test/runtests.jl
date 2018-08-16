# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

using BoundingSphere
const MB = BoundingSphere
using StaticArrays

if VERSION < v"0.7-"
    using Base.Test
    const stdout = STDOUT
    seed!(x) = srand(x)
else
    using Random
    using Random: seed!
    using LinearAlgebra
    using Test
    using InteractiveUtils: subtypes
end

include("helpers.jl")
include("test_broken.jl")
include("test_util.jl")
include("test_welzl.jl")
include("test_degenerate_examples.jl")
include("perf.jl")
