# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

@testset "2 points 2d" begin

    ball_ref = MB.SqBall{Array{Float64,1},Float64}([1.79979, -0.419288], 38.387992027461046)
    pts = StaticArrays.SArray{Tuple{2},Float64,1,2}[[0.379453, 0.65952], [1.29807, 1.06029], [1.66662, -0.225889]]

    c,r = boundingsphere(pts, WelzlPivot())
    ball = MB.SqBall(c, r^2)
    @test MB.allinside(pts, ball, rtol=1e-1)
    @test MB.allinside(pts, ball, rtol=1e-6)

    # not sure if WelzlMTF should give a good result
    c,r = boundingsphere(pts, WelzlMTF())
    ball = MB.SqBall(c, r^2)
    @test MB.allinside(pts, ball, rtol=1e-1)
    @test MB.allinside(pts, ball, rtol=1e-6)
end

@testset "2 points 1d" begin

    x = -5.
    y = x + 10*eps(Float64)
    pts = [[x], [y]]

    c,r = boundingsphere(pts, WelzlPivot())
    ball = MB.SqBall(c, r^2)
    @test MB.allinside(pts, ball, atol=100eps(Float64))

    # not sure if WelzlMTF should give a good result
    @test_broken boundingsphere(pts, WelzlMTF())
    # ball = MB.SqBall(c, r^2)
    # @test_broken MB.allinside(pts, ball, rtol=1e-1)
    # @test_broken MB.allinside(pts, ball, rtol=1e-6)
end

@testset "3 points 3d" begin


    ball_ref = MB.SqBall{Array{Float64,1},Float64}([1.05415, -2.52996, -0.584979], 20.953846606252085)
    pts = StaticArrays.SArray{Tuple{3},Float64,1,3}[[1.18024, 0.0853978, -3.01374], [0.20966, -3.69213, 3.76129], [4.51103, -0.881877, 1.92254]]

    c,r = boundingsphere(pts, WelzlPivot())
    ball = MB.SqBall(c, r^2)
    @test MB.allinside(pts, ball, rtol=1e-1)
    @test MB.allinside(pts, ball, rtol=1e-6)

    # not sure if WelzlMTF should give a good result
    c,r = boundingsphere(pts, WelzlMTF())
    ball = MB.SqBall(c, r^2)
    @test MB.allinside(pts, ball, rtol=1e-1)
    @test MB.allinside(pts, ball, rtol=1e-6)
end
