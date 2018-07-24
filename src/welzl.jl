# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

function welzl!(pts, bdry::BoundaryDevice, alg::WelzlMTF)
    bdry_len = length(bdry)

    support_count = 0

    ball = get_ball(bdry)
    if ismaxlength(bdry)
        support_count = 0
        return ball, support_count
    end

    for i in eachindex(pts)
        pt = pts[i]
        if !isinside(pt, ball)
            pts_i = prefix(pts, i-1)
            isstable = push_if_stable!(bdry, pt)
            if isstable
                ball, s = welzl!(pts_i, bdry, alg)
                @assert isinside(pt, ball, rtol=1e-2, atol=1e-10)
                pop!(bdry)
                move_to_front!(pts, i)
                support_count = s + 1
            end
        end
    end

    @assert bdry_len == length(bdry)
    ball, support_count
end

function find_max_excess(ball, pts, k1)
    T = eltype(first(pts))
    e_max = T(-Inf)
    k_max = k1 -1
    for k in k1:length(pts)
        pt = pts[k]
        e = sqdist(pt, center(ball)) - sqradius(ball)
        if  e > e_max
            e_max = e
            k_max = k
        end
    end
    e_max, k_max
end

function welzl!(pts, bdry, alg::WelzlPivot)
    t = 1
    alg_inner = WelzlMTF()
    @assert isempty(bdry)
    ball, s = welzl!(prefix(pts,t), bdry, alg_inner)
    for i in 1:alg.max_iterations
        @assert s <= t
        e, k = find_max_excess(ball, pts, t+1)

        P = eltype(pts)
        F = eltype(P)
        if e >  eps(F)  # TODO should this be a parameter of the algorithm?
            @assert t < k
            pt = pts[k]
            push_if_stable!(bdry, pt)
            ball_new, s_new = welzl!(prefix(pts,s), bdry, alg_inner)
            # @assert isinside(pt, ball_new, rtol=1e-6)
            pop!(bdry)
            @assert isempty(bdry)
            move_to_front!(pts,k)
            ball = ball_new
            t = s + 1
            s = s_new + 1
        else
            return ball
        end
    end
    return ball
end
