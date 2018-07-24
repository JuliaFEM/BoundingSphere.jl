# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

function prefix(pts, i)
    inds = eachindex(pts)
    index = first(inds) : i
    view(pts, index)
end

function move_to_front!(pts, i)
    @assert i in eachindex(pts)
    pt = pts[i]
    for j in eachindex(pts)
        qt = pts[j]
        pts[j] = pt
        pt = qt
        j == i && break
    end
    pts
end

function leq_approx(x,y;kw...)
    x < y || isapprox(x,y;kw...)
end
