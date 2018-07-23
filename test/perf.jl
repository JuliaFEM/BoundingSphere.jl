using BoundingSphere
const MB = BoundingSphere
using BenchmarkTools
using StaticArrays

for (dim, npoints) in [
                        (2, 100),
                        (3,1000)
                      ]
    pts = [@SVector(randn(3)) for _ in 1:npoints]
    for A in subtypes(MB.BoundingSphereAlg)
        alg = A()
        println("Running $alg on $npoints points of dimension $dim")
        trial = @benchmark boundingsphere($pts, $alg)
        show(STDOUT, MIME"text/plain"(), trial)
        println()
    end
end
