export miniball

export MiniballAlgorithm
export WelzlMTF
export WelzlPivot
export WelzlClassic
export Ritter


abstract type MiniballAlgorithm end

struct WelzlMTF <: MiniballAlgorithm end
struct WelzlClassic <: MiniballAlgorithm 
    shuffle::Bool
end
WelzlClassic(;shuffle=true) = WelzlClassic(shuffle)

struct WelzlPivot <: MiniballAlgorithm
    max_iterations::Int
end
function WelzlPivot(;max_iterations=1000)
    WelzlPivot(max_iterations)
end

needs_shuffle(alg) = false
needs_shuffle(alg::WelzlClassic) = alg.shuffle

function miniball!(pts, alg::WelzlMTF=WelzlMTF())
    bdry = create_boundary_device(pts, alg)
    ball, support_count = welzl!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

function miniball!(pts, alg::MiniballAlgorithm)
    if needs_shuffle(alg)
        shuffle!(pts)
    end
    bdry = create_boundary_device(pts, alg)
    ball = welzl!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

@noinline function miniball(pts, alg::MiniballAlgorithm=WelzlMTF())
    miniball!(copy(pts), alg)
end
