using MiniBallNext
const MB = MiniBallNext
using StaticArrays

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("test_util.jl")
include("test_welzl.jl")
