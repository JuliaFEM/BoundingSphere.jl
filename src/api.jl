export miniball

export WelzlMTF
export WelzlPivot
export Ritter


abstract type MiniballAlgorithm end

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
struct WelzlMTF <: MiniballAlgorithm end

"""
    WelzlPivot(;max_iterations=1000)

Welzl algorithm with pivoting. See Algorithm II in https://people.inf.ethz.ch/gaertner/subdir/texts/own_work/esa99_final.pdf.

## Pros
* Fast

## Cons
* In very rare cases can be numerically instable
"""
struct WelzlPivot <: MiniballAlgorithm
    max_iterations::Int
end

function WelzlPivot(;max_iterations=1000)
    WelzlPivot(max_iterations)
end

"""
    center, radius = miniball(pts [, algorithm=WelzlPivot()])

Compute the smallest sphere that contains each point in `pts`.

# Arguments

* pts: A list of points. Points should be vectors with floating point entries.
* algorithm: An optional algorithm to do the computation. See names(MiniBallNext) to get
"""
function miniball end

function miniball!(pts, alg::WelzlMTF=WelzlMTF())
    bdry = create_boundary_device(pts, alg)
    ball, support_count = welzl!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

function miniball!(pts, alg::MiniballAlgorithm)
    bdry = create_boundary_device(pts, alg)
    ball = welzl!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

@noinline function miniball(pts, alg::MiniballAlgorithm=WelzlMTF())
    miniball!(copy(pts), alg)
end
