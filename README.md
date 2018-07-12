# BoundingSphere

[![Build Status](https://travis-ci.org/JuliaFEM/BoundingSphere.jl.svg?branch=master)](https://travis-ci.org/JuliaFEM/BoundingSphere.jl)
[![codecov.io](https://codecov.io/github/JuliaFEM/BoundingSphere.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaFEM/BoundingSphere.jl?branch=master)
## Usage
```julia
using BoundingSphere

pts = [randn(3) for _ in 1:10]
center, radius = miniball(pts)

using StaticArrays
pts = [@SVector(randn(3)) for _ in 1:10] # use static arrays for performance
algorithm = Ritter() # fast but inaccurate
center, radius = miniball(pts, algorithm) # customize algorithm
```
