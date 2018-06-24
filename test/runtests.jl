using MiniBall
using StaticArrays

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

using MiniBall: SqBall, isinside, encloses

function create_ball_points(npoints, dim; p_boundary=1, p_rep=1/sqrt(npoints))
    P = SVector{dim, Float64}
    pts = P[]
    center = randn(dim)
    radius = 10*rand()
    ball = MiniBall.SqBall(center, radius^2)
    while length(pts) < npoints
        dir = randn(dim)
        r = if rand() < p_boundary
            radius
        else
            radius*rand()
        end
        pt = center + dir*r
        while true
            push!(pts, pt)
            (rand() >= p_rep) && break
            (length(pts) == npoints) && break
        end
    end
    shuffle!(pts)
    @assert length(pts) == npoints
    ball, pts
end

function random_test(npoints, dim; kw...)
    ball_ref, pts = create_ball_points(npoints, dim; kw...)
    c, r = miniball(pts)
    ball = SqBall(c, r^2)

    @test encloses(ball_ref, ball, rtol=1e-6)

    for pt in pts
        @test isinside(pt, ball, rtol=1e-6)
    end

end

@testset "random tests" begin
    srand(42)
    for dim in 1:10
        for npoints in 1:20
            random_test(npoints, dim, p_rep=0)
        end
    end
end
