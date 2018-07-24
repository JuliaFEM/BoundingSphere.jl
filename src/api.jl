# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

export boundingsphere

export WelzlMTF
export WelzlPivot
export Ritter


abstract type BoundingSphereAlg end

"""
    WelzlMTF()

Welzl algorithm with move to front heuristic.
See Algorithm I in https://people.inf.ethz.ch/gaertner/subdir/texts/own_work/esa99_final.pdf.
In almost all situations it is better to use [`WelzlPivot`](@ref) instead.

## Pros
* Fast for small examples

## Cons
* Prone to numerical stability issues
"""
struct WelzlMTF <: BoundingSphereAlg end

"""
    WelzlPivot(;max_iterations=1000)

Welzl algorithm with pivoting. See Algorithm II in https://people.inf.ethz.ch/gaertner/subdir/texts/own_work/esa99_final.pdf.

## Pros
* Fast

## Cons
* In very rare cases can be numerically instable
"""
struct WelzlPivot <: BoundingSphereAlg
    max_iterations::Int
end

function WelzlPivot(;max_iterations=1000)
    WelzlPivot(max_iterations)
end

"""
    center, radius = boundingsphere(pts [, algorithm=WelzlPivot()])

Compute the smallest sphere that contains each point in `pts`.

# Arguments

* pts: A list of points. Points should be vectors with floating point entries.
* algorithm: An optional algorithm to do the computation. See names(BoundingSphere) to get
"""
function boundingsphere end

function boundingsphere!(pts, alg::WelzlMTF=WelzlMTF())
    bdry = create_boundary_device(pts, alg)
    ball, support_count = welzl!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

function boundingsphere!(pts, alg::BoundingSphereAlg)
    bdry = create_boundary_device(pts, alg)
    ball = welzl!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

@noinline function boundingsphere(pts, alg::BoundingSphereAlg=WelzlMTF())
    boundingsphere!(copy(pts), alg)
end
