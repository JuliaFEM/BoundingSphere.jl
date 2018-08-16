# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/LICENSE

using StaticArrays

@testset "Generic type support" begin
    inputs = []
    pts_template = [[-1,0],[1,0],[0,-1],[0,1]]
    for F in subtypes(AbstractFloat)
        F == Float16 && continue  # Float16 fires assertions
        for V in [Vector{F}, SVector{2,F}]
            pts::Vector{V} = map(V, pts_template)
            push!(inputs, pts)
        end
    end
    for A in subtypes(MB.BoundingSphereAlg)
        for pts in inputs
            alg = A()
            c, r = boundingsphere(pts, alg)
            P = eltype(pts)
            F = eltype(P)
            @test typeof(c) == P
            @test typeof(r) == F
            @inferred boundingsphere(pts, alg)
            @test norm(c) < sqrt(eps(F))
            @test r ≈ 1
        end
    end
end

@testset "support_count" begin
    alg = WelzlMTF()

    a = @SVector [0.]
    b = @SVector [1.]
    c = @SVector [2.]

    pts = [a]
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 1

    pts = [a,b]
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 2

    pts = [a,b,c]
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 2
    @test Set(pts[1:support_count]) == Set([a,c])

    vals = randn(10)
    V = SVector{1, Float64}
    pts = map(V, vals)
    v1, v2 = map(V, extrema(vals))
    bdry = MB.create_boundary_device(pts, alg)
    ball, support_count = MB.welzl!(pts, bdry, alg)
    @test isempty(bdry)
    @test support_count == 2
    @test Set(pts[1:support_count]) == Set([v1,v2])

end

@testset "random support_count" begin
    seed!(42)
    alg = WelzlMTF()
    for dim in 1:10, npoints in 1:20
        ball, pts = create_ball_points(dim, npoints, p_rep=0, p_boundary=0)
        bdry = MB.create_boundary_device(pts, alg)
        ball, support_count = MB.welzl!(pts, bdry, alg)
        @test isempty(bdry)

        # support set suffices
        pts2 = pts[1:support_count]
        c, r = boundingsphere(pts2, alg)
        @test MB.center(ball) ≈ c
        @test MB.radius(ball) ≈ r

        # support set is not too large
        pts3 = copy(pts2)
        if length(pts3) > 1
            index = rand(1:length(pts3))
            deleteat!(pts3, index)
            c, r = boundingsphere(pts3, alg)
            @test r < MB.radius(ball)
        end
    end
end

@testset "random Ritter" begin
    seed!(42)
    for dim in 1:5
        for npoints in 1:100
            random_test(Ritter(), npoints, dim,
                        rtol_radius=0.3,
                        allow_broken=false
                       )
        end
    end
end

@testset "random WelzlPivot" begin

    @testset "random WelzlPivot small" begin
        seed!(42)
        alg = WelzlPivot()
        for dim in 1:3
            for npoints in 1:10
                    random_test(alg, npoints, dim,
                                allow_broken=false,
                                codim=0)
            end
        end
    end

    @testset "non degenerate" begin
        seed!(42)
        alg = WelzlPivot()
        for dim in 1:10
            for npoints in 1:100
                    random_test(alg, npoints, dim,
                                p_boundary=0,
                                p_rep=0,
                                codim=0,
                                allow_broken=false)
            end
        end
    end
    
    @testset "nasty" begin
        seed!(42)
        alg = WelzlPivot()
        for _ in 1:1
            for dim in 1:10
                for npoints in 1:100
                    random_test(alg, npoints, dim,
                                p_rep=1/sqrt(npoints),
                                p_boundary=0.5,
                                allow_broken=true)
    
                    random_test(alg, npoints, dim,
                                p_rep=1/sqrt(npoints),
                                p_boundary=1,
                                allow_broken=true)
    
                    random_test(alg, npoints, dim,
                                p_rep=0,
                                p_boundary=0,
                                allow_broken=true)
    
                    random_test(alg, npoints, dim,
                                allow_broken=true)

                    random_test(alg, npoints, dim,
                                codim=poisson(3),
                                allow_broken=true)
                end
            end
        end
    end
end

@testset "random WelzlMTF" begin
    seed!(42)
    for dim in 1:10, npoints in 1:20
        random_test(WelzlMTF(), npoints, dim,
                    codim=0, p_rep=0, p_boundary=0,
                    allow_broken=false)
    end
end
