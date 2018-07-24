# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

struct SqBall{P,F}
    center::P
    sqradius::F
end

function isinside(pt, ball::SqBall; atol=0, rtol=0)
    r2 = sqdist(pt, center(ball))
    R2 = sqradius(ball)
    r2 <= R2 || isapprox(r2, R2;atol=atol^2,rtol=rtol^2)
end

function allinside(pts, ball; kw...)
    for pt in pts
        isinside(pt, ball; kw...) || return false
    end
    true
end

center(b::SqBall) = b.center
radius(b::SqBall) = sqrt(b.sqradius)
sqradius(b::SqBall) = b.sqradius

dist(p1,p2) = norm(p1-p2)
sqdist(p1::AbstractVector, p2::AbstractVector) = sqnorm(p1-p2)
sqdist(x,y) = sqdist(y,x)

sqnorm(p) = sum(abs2,p)
