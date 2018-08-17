# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

bernoulli(p) = rand() <= p

function poisson(lambda)
    k = 0
    p = exp(-lambda)
    s = p
    x = rand()
    while x > s
        k += 1
        p = p * lambda / k
        s += p
    end
    k
end

struct Householder{T} <: AbstractMatrix{T}
    v::Vector{T}
end

function Base.size(h::Householder)
    l = length(h.v)
    l,l
end

function Base.getindex(h::Householder, i, j)
    I[i,j] - 2*h.v[i]*h.v[j]
end

function random_householder(dim)
    v = randn(dim)
    normalize!(v)
    Householder(v)
end

function random_orthogonal(dim)
    ret = if dim == 1
        [(-1.)^bernoulli(0.5)]
    else
        dim1 = dim - 1
        m = random_householder(dim)
        rot1 = random_orthogonal(dim1)
        m * [rot1 zeros(dim1, 1); zeros(1, dim1) 1]
    end
    @assert ret * ret' ≈ Matrix(I,dim,dim)
    ret
end

function random_embedding(dim_target, dim_src)
    inc = Matrix(I,dim_target, dim_src)
    rot = random_orthogonal(dim_target)
    M = SMatrix{dim_target, dim_src}
    M(rot * inc)
end

function create_ball_points(npoints, dim_src;
                            codim = poisson(0.3),
                            p_boundary=1,
                            p_rep=1/sqrt(npoints),
                            p_shuffle=0.5,
                           )
    dim_target = dim_src + codim
    @assert dim_src <= dim_target
    @assert 0 <= p_boundary <= 1
    @assert 0 <= p_rep <= 1
    @assert 0 <= p_shuffle <= 1
    F = Float64
    P = SVector{dim_src, F}
    pts = P[]
    center = randn(dim_src)
    R = 10*rand()
    ball = MB.SqBall(center, R^2)
    while length(pts) < npoints
        dir = normalize!(randn(dim_src))
        r = if bernoulli(p_boundary)
            R
        else
            R*rand()
        end
        pt = center + r*dir
        while true
            push!(pts, pt)
            bernoulli(p_rep) || break
            (length(pts) == npoints) && break
        end
    end
    bernoulli(p_shuffle) && shuffle!(pts)
    @assert length(pts) == npoints
    @assert eltype(pts) == P
    if dim_target > dim_src
        embedding = random_embedding(dim_target, dim_src)
        center = embedding * center
        pts = map(pt -> embedding*pt, pts)
        ball = MB.SqBall(center, R^2)
    end
    make_ball_containing_pts(ball, pts)
end

function make_ball_containing_pts(ball, pts)
    R2 = ball.sqradius
    center = ball.center
    for pt in pts
        R2pt = MB.sqdist(pt, center)
        if R2pt > R2
            @assert R2pt ≈ R2
            R2 = R2pt
        end
    end
    MB.SqBall(center, R2), pts
end

function random_test(alg, npoints, dim; 
                     rtol_inside=nothing,
                     atol_inside=nothing,
                     rtol_radius=nothing,
                     allow_broken::Bool=false,
                     kw...)
    ball_ref, pts = create_ball_points(npoints, dim; kw...)
    @assert MB.allinside(pts, ball_ref)
    P = eltype(pts)
    F = eltype(P)

    c, r = boundingsphere(pts, alg)
    ball = MB.SqBall(c, r^2)

    r_ref = MB.radius(ball_ref)
    if rtol_radius == nothing
        rtol_radius = sqrt(eps(max(r, r_ref)))
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

