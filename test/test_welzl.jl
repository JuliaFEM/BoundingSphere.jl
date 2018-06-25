using StaticArrays

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

function create_ball_points(npoints, dim; p_boundary=1, p_rep=1/sqrt(npoints))
    @assert 0 <= p_boundary <= 1
    @assert 0 <= p_rep <= 1
    F = Float64
    P = SVector{dim, F}
    pts = P[]
    center = randn(dim)
    R = 10*rand()
    ball = MiniBallNext.SqBall(center, R^2)
    while length(pts) < npoints
        dir = normalize!(randn(dim))
        r = if rand() < p_boundary
            R
        else
            R*rand()
        end
        @assert norm(dir) ≈ 1
        @assert r <= R
        pt = center + dir*r
        if !MB.isinside(pt, ball)
            R2new = MB.sqdist(center, pt)
            @assert R2new ≈ ball.sqradius
            @assert R2new >= ball.sqradius
            ball = MB.SqBall(center, R2new)
        end
        @assert MB.isinside(pt, ball)
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

function random_test(alg, npoints, dim; 
                     rtol_inside=nothing,
                     rtol_radius=nothing,
                     kw...)
    ball_ref, pts = create_ball_points(npoints, dim; kw...)
    c, r = miniball(pts, alg)
    ball = MB.SqBall(c, r^2)

    r_ref = MB.radius(ball_ref)
    if rtol_radius == nothing
        rtol_radius = Base.rtoldefault(r, r_ref)
    end
    @test r <= r_ref || isapprox(r, r_ref; rtol=rtol_radius)

    @assert MB.allinside(pts, ball_ref)
    if rtol_inside == nothing
        rtol_inside = 1e-6
    end
    @test MB.allinside(pts, ball, rtol=rtol_inside)

end

@testset "random Ritter" begin
    srand(42)
    for dim in 1:5
        for npoints in 1:100
            random_test(Ritter(), npoints, dim,
                        rtol_radius=0.2)
        end
    end
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

