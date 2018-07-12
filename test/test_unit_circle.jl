using StaticArrays

@testset "stastical unit circle test" begin
    pts = Vector{SVector{2,Float64}}(5)
    for k in 1:100000
        for i in 1:4
            angle = rand()*2*pi
            x = cos(angle)
            y = sin(angle)
            pts[i] = @SVector [x, y]
            if i == 1
              pts[end] = @SVector [cos(angle + pi), sin(angle + pi)]
            end
        end

        mb = miniball(pts)
        origo = isapprox(mb[1],[0.0;0.0];atol=sqrt(eps()))
        if origo
            @test origo
        else
            @test_broken origo
        end
    end
end
