__precompile__()
module MiniBallNext

using ArgCheck # TODO get rid of this

include("api.jl")
include("geometry.jl")
include("welzl.jl")
include("util.jl")

end # module
