# BoundingSphere

[![][gitter-img]][gitter-url]
[![][travis-img]][travis-url]
[![][pkg-0.6-img]][pkg-0.6-url]
[![][pkg-0.7-img]][pkg-0.7-url]
[![][coveralls-img]][coveralls-url]
[![][docs-stable-img]][docs-stable-url]
[![][docs-latest-img]][docs-latest-url]
[![][issues-img]][issues-url]
[![][appveyor-img]][appveyor-url]

Package contains algorithms to calculate smallest enclosing sphere for a given
set of points in N dimensions.

BoundingSphere.jl is a complete rewrite from scratch of [Miniball.jl](https://github.com/JuliaFEM/Miniball.jl).
See [Miniball.jl issue #28](https://github.com/JuliaFEM/Miniball.jl/issues/28).

## Usage

```julia
using BoundingSphere

pts = [randn(3) for _ in 1:10]
center, radius = boundingsphere(pts)

using StaticArrays
pts = [@SVector(randn(3)) for _ in 1:10] # use static arrays for performance
algorithm = Ritter() # fast but inaccurate
center, radius = boundingsphere(pts, algorithm) # customize algorithm
```
[gitter-img]: https://badges.gitter.im/Join%20Chat.svg
[gitter-url]: https://gitter.im/JuliaFEM/JuliaFEM.jl

[contrib-url]: https://juliafem.github.io/BoundingSphere.jl/latest/man/contributing/
[discourse-tag-url]: https://discourse.julialang.org/tags/boundingsphere

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://juliafem.github.io/BoundingSphere.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliafem.github.io/BoundingSphere.jl/stable

[travis-img]: https://travis-ci.org/JuliaFEM/BoundingSphere.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaFEM/BoundingSphere.jl

[coveralls-img]: https://coveralls.io/repos/github/JuliaFEM/BoundingSphere.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaFEM/BoundingSphere.jl?branch=master

[issues-img]: https://img.shields.io/github/issues/JuliaFEM/BoundingSphere.jl.svg
[issues-url]: https://github.com/JuliaFEM/BoundingSphere.jl/issues

[pkg-0.6-img]: http://pkg.julialang.org/badges/BoundingSphere_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=BoundingSphere&ver=0.6
[pkg-0.7-img]: http://pkg.julialang.org/badges/BoundingSphere_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=BoundingSphere&ver=0.7

[appveyor-img]: https://ci.appveyor.com/api/projects/status/s1vk9v0sxbmr2pen/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaFEM/boundingsphere-jl/branch/master
