struct WelzlMTF <: MiniballAlgorithm end

const F = Float64
const P = Vector{F}

struct ProjectorStack 
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

abstract type BoundaryDevice end
mutable struct GaertnerBdry <: BoundaryDevice
    centers::Vector{P}
    square_radii::Vector{F} # square radii
    projector::ProjectorStack
end

function create_boundary_device(pts, alg::WelzlMTF)
    projector = ProjectorStack([])
    GaertnerBdry([],[],projector)
end

function npoints(b::GaertnerBdry)
    @assert length(b.centers) ==
    length(b.square_radii)

    return length(b.centers)
end

function Base.push!(b::GaertnerBdry, pt)
    if npoints(b) == 0
        push!(b.square_radii, zero(eltype(pt)))
        push!(b.centers, pt)
        dim = length(pt)
        return b
    end

    q0 = first(b.centers)
    center = b.centers[end]
    C  = center - q0
    r2 = b.square_radii[end]

    Qm = pt - q0
    M = b.projector
    Qm_bar = M*Qm

    e = sqdist(Qm, C) - r2
    z = 2*sqdist(Qm, Qm_bar)

    @assert abs(z) > eps(F) # TODO small z

    center_new  = center + (e/z) * (Qm - Qm_bar)
    r2new = r2 + (e^2)/(2z)

    residue = Qm - Qm_bar
    residue_norm = residue / norm(residue)
    push!(b.projector, residue / norm(residue))

    push!(b.centers, center_new)
    push!(b.square_radii, r2new)
    
    b
end

function Base.pop!(b::GaertnerBdry)
    @argcheck npoints(b) > 0
    pop!(b.centers)
    pop!(b.square_radii)
    if npoints(b) > 0 # npoints has been poped here
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

function prefix(pts, i)
    inds = eachindex(pts)
    index = first(inds) : i
    view(pts, index)
end

function move_to_front!(pts, i)
    pt = pts[i]
    for j in prefix(pts, i)
        qt = pts[i]
        pts[i] = pt
        pt = qt
    end
    pts
end

function mb!(pts, bdry::BoundaryDevice, alg::WelzlMTF)

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
    if npoints(bdry) == dim + 1
        return ball
    end

    for i in eachindex(pts)
        pt = pts[i]
        if !isinside(pt, ball)
            pts_i = prefix(pts, i-1)
            push!(bdry, pt)
            ball = mb!(pts_i, bdry, alg)
            @assert isinside(pt, ball, rtol=1e-2)
            pop!(bdry)
            move_to_front!(pts, i)
        end
    end
    ball
end
