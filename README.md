# MiniBallNext

[![Build Status](https://travis-ci.org/jw3126/MiniBallNext.jl.svg?branch=master)](https://travis-ci.org/jw3126/MiniBallNext.jl)
[![codecov.io](https://codecov.io/github/jw3126/MiniBallNext.jl/coverage.svg?branch=master)](http://codecov.io/github/jw3126/MiniBallNext.jl?branch=master)
## Usage
```julia
using MiniBallNext

pts = [randn(3) for _ in 1:10]
center, radius = miniball(pts)

using StaticArrays
pts = [@SVector(randn(3)) for _ in 1:10] # use static arrays for performance
algorithm = Ritter() # fast but inaccurate
center, radius = miniball(pts, algorithm) # customize algorithm
```
