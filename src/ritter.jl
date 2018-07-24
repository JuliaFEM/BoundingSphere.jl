# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

"""
    Ritter()

## Pros
* extremly fast
* simple

## Cons
* Very inaccurate.
"""
struct Ritter <: BoundingSphereAlg end

function max_distance_point(pts, pt1)
    pt_best = first(pts)
    d2_best = sqdist(pt_best, pt1)
    for pt in pts
        d2 = sqdist(pt, pt1)
        if d2 > d2_best
            pt_best = pt
            d2_best = d2
        end
    end
    pt_best
end

function sqsphere_two_points(pt1,pt2)
    c = map(middle, pt1, pt2)
    r2 = sqdist(pt1, pt2) / 4
    c, r2
end

@noinline function ritter(pts)
    pt1 = first(pts)
    pt2 = max_distance_point(pts, pt1)
    c, r2 = sqsphere_two_points(pt1, pt2)
    for pt in pts
        if sqdist(pt, c) > r2
            direction = (c - pt) / norm(c - pt)
            r = sqrt(r2)
            pt_op = c + direction * r
            c,r2 = sqsphere_two_points(pt, pt_op)
        end
    end
    r = sqrt(r2)
    c, r
end

boundingsphere(pts, alg::Ritter) = ritter(pts)
