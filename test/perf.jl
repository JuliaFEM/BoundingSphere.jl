using MiniBallNext
using MiniBallNext
const MB = MiniBallNext
using BenchmarkTools
using StaticArrays

for (dim, npoints) in [
                        (2, 100),
                        (3,1000)
                      ]
    pts = [@SVector(randn(3)) for _ in 1:npoints]
    for A in subtypes(MB.MiniballAlgorithm)
        alg = A()
        println("Running $alg on $npoints points of dimension $dim")
        trial = @benchmark miniball($pts, $alg)
        show(STDOUT, MIME"text/plain"(), trial)
        println()
    end
end
