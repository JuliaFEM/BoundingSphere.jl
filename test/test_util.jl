# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

@testset "prefix" begin
    @test MB.prefix(["a", "b", "c", "d"], 0) == []
    @test MB.prefix(["a", "b", "c", "d"], 1) == ["a"]
    @test MB.prefix(["a", "b", "c", "d"], 2) == ["a","b"]

    # mut
    data = [1,2,3]
    p = MB.prefix(data, 1)
    p[1] = 42
    @test data == [42,2,3]
end

@testset "move_to_front!" begin
    pts_initial = ["a", "b", "c", "d"]

    pts = deepcopy(pts_initial)
    MB.move_to_front!(pts, 1)
    @test pts == pts_initial

    pts = deepcopy(pts_initial)
    MB.move_to_front!(pts, 2)
    @test pts == ["b", "a", "c", "d"]

    pts = deepcopy(pts_initial)
    MB.move_to_front!(pts, 4)
    @test pts == ["d", "a", "b", "c"]
end

@testset "isinside" begin
    @test MB.isinside([0.], MB.SqBall([0.], 0.))
    @test !MB.isinside([0.], MB.SqBall([NaN], 0.))
    @test !MB.isinside([0.], MB.SqBall([0.], NaN))

end
