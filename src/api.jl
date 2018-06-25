export miniball

export MiniballAlgorithm
export WelzlMTF
export WelzlPivot
export Ritter


abstract type MiniballAlgorithm end

struct WelzlMTF <: MiniballAlgorithm end

struct WelzlPivot <: MiniballAlgorithm
    max_iterations::Int
end
function WelzlPivot(;max_iterations=1000)
    WelzlPivot(max_iterations)
end

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
