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
        # pts = map(pt -> map(BigFloat,pt), pts)
        # @show length(pts)
        # @show length(first(pts))
        # @show ball_ref
        # @show pts
        # @show ball
    else
        @test contains_all_points
    end
end

