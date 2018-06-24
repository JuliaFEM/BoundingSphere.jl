export miniball
export WelzlMTF
export WelzlPivot

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
    ball, support_count = mb!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

function miniball!(pts, alg::MiniballAlgorithm)
    ball = mb!(pts, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

function miniball(pts, alg::MiniballAlgorithm=WelzlMTF())
    miniball!(copy(pts), alg)
end
