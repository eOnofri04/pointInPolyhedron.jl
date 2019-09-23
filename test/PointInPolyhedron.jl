using Test

@testset "Determine Simplex Solid Angle" begin
    p = [0.0; 0.0; 0.0];
    x = [1.0; 0.0; 0.0];
    y = [0.0; 1.0; 0.0];
    z = [0.0; 0.0; 1.0];

    @testset "Octants" begin

        @test isapprox(simplexSolidAngle(p, [+x +y +z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [+x +y -z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [+x -y +z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [+x -y -z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [-x +y +z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [-x +y -z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [-x -y +z])/π, 0.5);
        @test isapprox(simplexSolidAngle(p, [-x -y -z])/π, 0.5);

    end

    @testset "Hexant" begin
        @testset "xy-Hexant" begin

            @test isapprox(simplexSolidAngle(p, [+x +x+y +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x+y +y +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+y +y-x +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+y-x -x +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x -x-y +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x-y -y +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-y -y+x +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-y+x +x +z])/π, 0.25);

            @test isapprox(simplexSolidAngle(p, [+x +x+y -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x+y +y -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+y +y-x -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+y-x -x -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x -x-y -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x-y -y -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-y -y+x -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-y+x +x -z])/π, 0.25);
        end

        @testset "xz-Hexant" begin
            @test isapprox(simplexSolidAngle(p, [+x +y +x+z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x+z +y +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+z +y +z-x])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+z-x +y -x])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x +y -x-z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x-z +y -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-z +y -z+x])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-z+x +y +x])/π, 0.25);

            @test isapprox(simplexSolidAngle(p, [+x -y +x+z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x+z -y +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+z -y +z-x])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+z-x -y -x])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x -y -x-z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-x-z -y -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-z -y -z+x])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [-z+x -y +x])/π, 0.25);
        end

        @testset "yz-Hexant" begin
            @test isapprox(simplexSolidAngle(p, [+x +y +y+z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x +y+z +z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x +z +z-y])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x +z-y -y])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x -y -y-z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x -y-z -z])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x -z -z+y])/π, 0.25);
            @test isapprox(simplexSolidAngle(p, [+x -z+y +y])/π, 0.25);
        end
    end
end

@testset "Point of a Cube" begin
    V = [
        0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
        0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
        0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0
    ]
    EV = [
        [1, 2], [2, 4], [4, 3], [3, 1],
        [5, 6], [6, 8], [8, 7], [7, 5],
        [1, 5], [2, 6], [4, 8], [3, 5]
    ];
    FV = [
        [1, 3, 4, 2], [5, 6, 8, 7],
        [1, 5, 7, 3], [2, 4, 8, 6],
        [1, 2, 6, 5], [3, 7, 8, 4],
    ];
    triFV = [
        [1, 3, 4], [4, 2, 1], [5, 6, 8], [8, 7, 5],
        [1, 5, 7], [7, 3, 1], [2, 4, 8], [8, 6, 2],
        [1, 2, 6], [6, 5, 1], [3, 7, 8], [8, 4, 3],
    ];

    @testset "Internal Points" begin
        @test pointInPolyhedron3D([0.5; 0.5; 0.5], V, EV, triFV) == 1;
        @test pointInPolyhedron3D([0.01; 0.01; 0.01], V, EV, triFV) == 1;
        @test pointInPolyhedron3D([0.01; 0.01; 0.00000001], V, EV, triFV) == 1;
        @test pointInPolyhedron3D([0.9999999; 0.5; 10^-5], V, EV, triFV) == 1;
    end

    @testset "Outer Points" begin
        @test pointInPolyhedron3D([0.01; 0.01; -0.0000001], V, EV, triFV) == -1;
        @test pointInPolyhedron3D([0.5; 0.5; 1.5], V, EV, triFV) == -1;
        @test pointInPolyhedron3D([0.5; 0.5; 1.0000000001], V, EV, triFV) == -1;
        @test pointInPolyhedron3D([0.0; 0.0; -0.00000001], V, EV, triFV) == -1;
    end

    @testset "Edge Points" begin
        @test pointInPolyhedron3D([0.5; 0.5; 0.0], V, EV, triFV) == 0;
        @test pointInPolyhedron3D([0.5; 0.0; 0.0], V, EV, triFV) == 0;
        @test pointInPolyhedron3D([0.0; 0.0; 0.0], V, EV, triFV) == 0;
        @test pointInPolyhedron3D([1.0; 0.5; 0.5], V, EV, triFV) == 0;
    end
end
