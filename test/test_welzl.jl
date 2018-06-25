using StaticArrays

function create_ball_points(npoints, dim; p_boundary=1, p_rep=1/sqrt(npoints))
    @assert 0 <= p_boundary <= 1
    @assert 0 <= p_rep <= 1
    P = SVector{dim, Float64}
    pts = P[]
    center = randn(dim)
    radius = 10*rand()
    ball = MiniBallNext.SqBall(center, radius^2)
    while length(pts) < npoints
        dir = normalize!(randn(dim))
        r = if rand() < p_boundary
            radius
        else
            radius*rand()
        end
        pt = center + dir*r
        @assert MB.isinside(pt, ball, rtol=1e-6)
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

function random_test(alg, npoints, dim; rtol=1e-6, kw...)
    ball_ref, pts = create_ball_points(npoints, dim; kw...)
    c, r = miniball(pts, alg)
    ball = MB.SqBall(c, r^2)

    r_ref = MB.radius(ball_ref)
    @test r <= r_ref || r ≈ r_ref

    @assert MB.allinside(pts, ball_ref, rtol=rtol)
    @test MB.allinside(pts, ball, rtol=rtol)

end

@testset "random WelzlMTF" begin
    srand(42)
    for dim in 1:10, npoints in 1:20
        random_test(WelzlMTF(), npoints, dim, p_rep=0, p_boundary=0)
    end
end

@testset "random WelzlPivot" begin
    srand(42)
    for dim in 1:10
        for npoints in 1:100
            random_test(WelzlPivot(), npoints, dim,
                        p_rep=1/sqrt(npoints),
                        p_boundary=0.5)
        end
    end
end

@testset "Generic type support" begin
    inputs = []
    for F in subtypes(AbstractFloat)
        F == Float16 && continue  # Float16 fires assertions
        pts = Vector{F}[F.(randn(3)) for _ in 1:5]
        push!(inputs, pts)
        pts = [F.(@SVector(randn(3))) for _ in 1:5]
        push!(inputs, pts)
    end
    for A in subtypes(MB.MiniballAlgorithm)
        for pts in inputs
            alg = A()
            c, r = miniball(pts, alg)
            P = eltype(pts)
            F = eltype(P)
            @test typeof(c) == P
            @test typeof(r) == F
            @inferred miniball(pts, alg)
        end
    end
end
