struct SqBall{C,R}
    center::C
    sqradius::R
end

function isinside(pt, ball::SqBall; atol=0, rtol=0)
    r2 = sqdist(pt, ball)
    R2 = sqradius(ball)
    r2 <= R2 || isapprox(r2, R2;atol=atol,rtol=rtol)
end

center(b::SqBall) = b.center
radius(b::SqBall) = sqrt(b.sqradius)
sqradius(b::SqBall) = b.sqradius

dist(p1,p2) = norm(p1-p2)
sqdist(p1::AbstractVector, p2::AbstractVector) = sqnorm(p1-p2)
sqdist(x,y) = sqdist(y,x)
sqdist(pt,b::SqBall) = sqdist(pt, b.center)

sqnorm(p) = sum(abs2,p)

function encloses(bigger::SqBall, smaller::SqBall;atol=0,rtol=0)
    r = dist(center(bigger), center(smaller)) + radius(smaller)
    R = radius(bigger)
    r <= R || isapprox(r,R;atol=atol,rtol=rtol)
end
