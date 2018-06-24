export miniball

abstract type MiniballAlgorithm end

function miniball!(pts, alg::MiniballAlgorithm=WelzlMTF())
    bdry = create_boundary_device(pts, alg)
    ball = mb!(pts, bdry, alg)
    r = radius(ball)
    c = center(ball)
    c, r
end

function miniball(pts, alg::MiniballAlgorithm=WelzlMTF())
    miniball!(copy(pts), alg)
end
