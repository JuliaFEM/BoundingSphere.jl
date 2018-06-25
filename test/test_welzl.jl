using StaticArrays

@testset "Generic type support" begin
    inputs = []
    pts_template = [[-1,0],[1,0],[0,-1],[0,1]]
    for F in subtypes(AbstractFloat)
        F == Float16 && continue  # Float16 fires assertions
        for V in [Vector{F}, SVector{2,F}]
            pts::Vector{V} = map(V, pts_template)
            push!(inputs, pts)
        end
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
            @test norm(c) < sqrt(eps(F))
            @test r ≈ 1
        end
    end
end

bernoulli(p) = rand() <= p

function create_ball_points(npoints, dim; 
                            p_boundary=1, 
                            p_rep=1/sqrt(npoints),
                            p_shuffle=0.5,
                           )
    @assert 0 <= p_boundary <= 1
    @assert 0 <= p_rep <= 1
    @assert 0 <= p_shuffle <= 1
    F = Float64
    P = SVector{dim, F}
    pts = P[]
    center = randn(dim)
    R = 10*rand()
    ball = MiniBallNext.SqBall(center, R^2)
    while length(pts) < npoints
        dir = normalize!(randn(dim))
        r = if bernoulli(p_boundary)
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
            bernoulli(p_rep) || break
            (length(pts) == npoints) && break
        end
    end
    bernoulli(p_shuffle) && shuffle!(pts)
    @assert length(pts) == npoints
    @assert eltype(pts) == P
    ball, pts
end

function random_test(alg, npoints, dim; 
                     rtol_inside=nothing,
                     atol_inside=nothing,
                     rtol_radius=nothing,
                     allow_broken::Bool=true,
                     kw...)
    ball_ref, pts = create_ball_points(npoints, dim; kw...)
    @assert MB.allinside(pts, ball_ref)
    P = eltype(pts)
    F = eltype(P)

    c, r = miniball(pts, alg)
    ball = MB.SqBall(c, r^2)

    r_ref = MB.radius(ball_ref)
    if rtol_radius == nothing
        rtol_radius = Base.rtoldefault(r, r_ref)
    end
    issmall = r <= r_ref || isapprox(r, r_ref; rtol=rtol_radius)
    if allow_broken && !issmall
        @test_broken issmall
    else
        @test issmall
    end

    if rtol_inside == nothing
        rtol_inside = 1e-6
    end
    if atol_inside == nothing
        atol_inside = 100eps(F)
    end
    contains_all_points = MB.allinside(pts, ball, rtol=rtol_inside, atol=atol_inside)
    if allow_broken && !contains_all_points
        @test_broken contains_all_points
        @show length(pts)
        @show length(first(pts))
        @show ball_ref
        @show pts
        @show ball
    else
        @test contains_all_points
    end
end

@testset "support_count" begin
    alg = WelzlMTF()

    a = @SVector [0.]
    b = @SVector [1.]
    c = @SVector [2.]

    pts = [a]
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 1

    pts = [a,b]
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 2

    pts = [a,b,c]
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 2
    @test Set(pts[1:support_count]) == Set([a,c])

    vals = randn(10)
    V = SVector{1, Float64}
    pts = map(V, vals)
    v1, v2 = map(V, extrema(vals))
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 2
    @test Set(pts[1:support_count]) == Set([v1,v2])

end

@testset "random support_count" begin
    srand(42)
    alg = WelzlMTF()
    for dim in 1:10, npoints in 1:20
        ball, pts = create_ball_points(dim, npoints, p_rep=0, p_boundary=0)
        bdry = MB.create_boundary_device(pts, alg)
        ball, support_count = MB.welzl!(pts, bdry, alg)
        @test isempty(bdry)

        # support set suffices
        pts2 = pts[1:support_count]
        c, r = miniball(pts2, alg)
        @test MB.center(ball) ≈ c
        @test MB.radius(ball) ≈ r

        # support set is not too large
        pts3 = copy(pts2)
        if length(pts3) > 1
            index = rand(1:length(pts3))
            deleteat!(pts3, index)
            c, r = miniball(pts3, alg)
            @test r < MB.radius(ball)
        end
    end
end

@testset "random Ritter" begin
    srand(42)
    for dim in 1:5
        for npoints in 1:100
            random_test(Ritter(), npoints, dim,
                        rtol_radius=0.3,
                       )
        end
    end
end

@testset "random WelzlPivot" begin
    srand(42)
    for _ in 1:1
        for dim in 1:10
            for npoints in 1:100
                random_test(WelzlPivot(), npoints, dim,
                            p_rep=1/sqrt(npoints),
                            p_boundary=0.5)

                random_test(WelzlPivot(), npoints, dim,
                            p_rep=1/sqrt(npoints),
                            p_boundary=1)

                random_test(WelzlPivot(), npoints, dim,
                            p_rep=0,
                            p_boundary=0)

                random_test(WelzlPivot(), npoints, dim,
                           )
            end
        end
    end
end

@testset "random WelzlMTF" begin
    srand(42)
    for dim in 1:10, npoints in 1:20
        random_test(WelzlMTF(), npoints, dim, p_rep=0, p_boundary=0)
    end
end

@testset "random WelzlClassic" begin
    srand(42)
    alg = WelzlClassic()
    for dim in 1:10, npoints in 1:20
        random_test(alg, npoints, dim, p_rep=0, p_boundary=0)
    end
end

