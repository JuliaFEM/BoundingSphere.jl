struct ProjectorStack{P <: AbstractVector}
    # matrix that is decomposed into Σ v_i ⊗ v_i* for
    # an orthonormal system v_i
    vs::Vector{P}
end

function Base.push!(p::ProjectorStack, v)
    @assert norm(v) ≈ 1
    push!(p.vs, v)
    p
end
function Base.pop!(p::ProjectorStack) 
    pop!(p.vs)
    p
end

function Base.:*(p::ProjectorStack, v::AbstractVector)
    ret = zero(v)
    for vi in p.vs
        ret = ret + vi*dot(vi, v)
    end
    ret
end

"""
    BoundaryDevice

Finds unique spheres determined by prescribed affine independent
boundary points. In the welzl algorithm this problem needs to be solved
in series, where points are pushed and popped from to the boundary.

Subtypes must implement the following interface:

* push_if_stable!(device, pt)::Bool :
* pop!(device): Remove last point from the boundary.
* get_ball(device)::SqBall : Get the last ball from the device.
"""
abstract type BoundaryDevice end

"""
    GaertnerBdry

BoundaryDevice that corresponds to M_B in Section 4 of Gaertners paper.

See also: [BoundaryDevice](@ref)
"""
mutable struct GaertnerBdry{P<:AbstractVector,
                            F<:AbstractFloat} <: BoundaryDevice
    centers::Vector{P}
    square_radii::Vector{F} # square radii
    projector::ProjectorStack{P}
end

function create_boundary_device(pts, alg)
    P = eltype(pts)
    F = eltype(P)
    projector = ProjectorStack(P[])
    centers = P[]
    square_radii = F[]
    GaertnerBdry(centers, square_radii, projector)
end

function npoints(b::GaertnerBdry)
    @assert length(b.centers) ==
    length(b.square_radii)

    return length(b.centers)
end

function push_if_stable!(b::GaertnerBdry, pt)
    if npoints(b) == 0
        push!(b.square_radii, zero(eltype(pt)))
        push!(b.centers, pt)
        dim = length(pt)
        return true
    end

    q0 = first(b.centers)
    center = b.centers[end]
    C  = center - q0
    r2 = b.square_radii[end]

    Qm = pt - q0
    M = b.projector
    Qm_bar = M*Qm
    residue = Qm - Qm_bar

    e = sqdist(Qm, C) - r2
    z = 2*sqnorm(residue)

    # should we use norm(residue) instead of z here?
    # seems more intuitive, OTOH z is used in the paper
    isstable = abs(z) > eps(eltype(pt))
    if isstable

        center_new  = center + (e/z) * residue
        r2new = r2 + (e^2)/(2z)

        push!(b.projector, residue / norm(residue))
        push!(b.centers, center_new)
        push!(b.square_radii, r2new)
    end
    isstable
end

function Base.pop!(b::GaertnerBdry)
    n = npoints(b)
    pop!(b.centers)
    pop!(b.square_radii)
    if n >= 2
        pop!(b.projector)
    end
    b
end

function get_ball(b::GaertnerBdry)
    @argcheck npoints(b) > 0
    c = b.centers[end]
    r2 = b.square_radii[end]
    SqBall(c,r2)
end

function welzl!(pts, bdry::BoundaryDevice, alg::WelzlMTF)

    dim = if isempty(pts)
        length(length(first(bdry.centers)))
    else
        length(first(pts))
    end

    if npoints(bdry) == 0
        c = first(pts)
        r2 = zero(eltype(c))
        ball = SqBall(c,r2)
    else
        ball = get_ball(bdry)
    end

    support_count = npoints(bdry)
    if support_count == dim + 1
        return ball, support_count
    end

    for i in eachindex(pts)
        pt = pts[i]
        if !isinside(pt, ball)
            pts_i = prefix(pts, i-1)
            isstable = push_if_stable!(bdry, pt)
            if isstable
                ball, support_count = welzl!(pts_i, bdry, alg)
                @assert isinside(pt, ball, rtol=1e-2)
                pop!(bdry)
                move_to_front!(pts, i)
            end
        end
    end
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

function welzl!(pts, alg::WelzlPivot)
    t = 1
    alg_inner = WelzlMTF()
    bdry = create_boundary_device(pts, alg)
    @assert npoints(bdry) == 0
    ball, s = welzl!(prefix(pts,t), bdry, alg_inner)
    for i in 1:alg.max_iterations
        e, k = find_max_excess(ball, pts, t+1)
        # @assert s <= t
        if e > 0
            @assert t < k
            pt = pts[k]
            push_if_stable!(bdry, pt)
            ball, s2 = welzl!(prefix(pts,s), bdry, alg_inner)
            @assert isinside(pt, ball, rtol=1e-2)

            pop!(bdry)
            @assert npoints(bdry) == 0
            move_to_front!(pts,k)
            t = s + 1
            s = s2 + 1
        else
            return ball
        end
    end
    return ball
end
