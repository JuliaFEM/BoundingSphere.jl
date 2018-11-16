# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

using Documenter
using BoundingSphere

makedocs(modules=[BoundingSphere],
         format = :html,
         checkdocs = :all,
         sitename = "BoundingSphere.jl",
         analytics = "UA-83590644-1",
         pages = ["index.md", "api.md"])
