# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

__precompile__()
module BoundingSphere

include("api.jl")
include("boundary.jl")
include("geometry.jl")
include("welzl.jl")
include("ritter.jl")
include("util.jl")

end # module
