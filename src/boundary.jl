# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

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
* length(device)::Int : Get the current count of boundary points
* ismaxlength(device)::Bool: Check if there are dim+1 boundary points in the device
"""
abstract type BoundaryDevice end

Base.isempty(b::BoundaryDevice) = length(b) == 0

"""
    GaertnerBdry

BoundaryDevice that corresponds to M_B in Section 4 of Gaertners paper.

See also: [BoundaryDevice](@ref)
"""
mutable struct GaertnerBdry{P<:AbstractVector,
                            F<:AbstractFloat} <: BoundaryDevice
    centers::Vector{P}
    square_radii::Vector{F}
    # projection onto of affine space spanned by points
    # shifted such that first point becomes origin
    projector::ProjectorStack{P}
    empty_center::P # center of ball spanned by empty boundary
end

function create_boundary_device(pts, alg)
    P = eltype(pts)
    F = eltype(P)
    projector = ProjectorStack(P[])
    centers = P[]
    square_radii = F[]
    empty_center = F(NaN)*first(pts)
    GaertnerBdry(centers, square_radii, projector, empty_center)
end

function Base.length(b::GaertnerBdry)
    @assert length(b.centers) ==
    length(b.square_radii)

    return length(b.centers)
end

function push_if_stable!(b::GaertnerBdry, pt)
    if isempty(b)
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

    # TODO should we use norm(residue) instead of z here?
    # seems more intuitive, OTOH z is used in the paper
    tol = eps(eltype(pt)) * max(r2, one(r2))
    @assert !isnan(tol)
    isstable = abs(z) > tol
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
    n = length(b)
    pop!(b.centers)
    pop!(b.square_radii)
    if n >= 2
        pop!(b.projector)
    end
    b
end

function get_ball(b::GaertnerBdry)
    if isempty(b)
        c = b.empty_center
        r2 = zero(eltype(c))
    else
        c = b.centers[end]
        r2 = b.square_radii[end]
    end
    SqBall(c,r2)
end

function ismaxlength(b::GaertnerBdry)
    dim = length(b.empty_center)
    length(b) == dim + 1
end
