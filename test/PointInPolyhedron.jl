using Test;
include("../src/pointInPolyhedron.jl");


@testset "Other Geometric Stuff" begin
    @testset "areCoplanar" begin
        x = [1.0; 0.0; 0.0];
        y = [0.0; 1.0; 0.0];
        z = [0.0; 0.0; 1.0];

        @test areCoplanar(x, x, y, z) == true;
        @test areCoplanar(x, y, y, z) == true;
        @test areCoplanar(x, y, z, z) == true;
        @test areCoplanar(x, 2*x, y, 2*y) == true;
        @test areCoplanar(x, 2*x, 3*x, y) == true;
        @test areCoplanar(x, 2*x, y, x+y) == true;
        @test areCoplanar(x, x+y, y, 2*x+3*y) == true;
        @test areCoplanar(x, z, x+z, 2*x-z) == true;

        @test areCoplanar(x, y, z, x+2*y) == false;
        @test areCoplanar(x, y, z, x-x) == false;
        @test areCoplanar(x, y, x+y, z) == false;
        @test areCoplanar(x, y, x+y, z*10e-7) == false;
        @test areCoplanar(x, y, x-y, z*10e-8) == false;
        @test areCoplanar(x, y, z, (1-10e-8)*z) == false;
        @test areCoplanar(x, y, 2*x, x-z) == false;
        @test areCoplanar(x, 2*x, 2*y, z) == false;
    end
end


@testset "Determine Simplex Solid Angle" begin
    p = [0.0; 0.0; 0.0];
    x = [1.0; 0.0; 0.0];
    y = [0.0; 1.0; 0.0];
    z = [0.0; 0.0; 1.0];

    @testset "Quadrants" begin
        @test isapprox(simplexSolidAngle(p, [+x -x +z])/π, 1.0);
        @test isapprox(simplexSolidAngle(p, [+x -x+y*1e-10 +z])/π, 1.0);
    end

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

    @testset "Other Angles" begin
        @test isapprox(simplexSolidAngle(p, [+x +y-x +z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [+x +y-x -z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [+x -y-x +z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [+x -y-x -z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [-x +y+x +z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [-x +y+x -z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [-x -y+x +z])/π, 0.75)
        @test isapprox(simplexSolidAngle(p, [-x -y+x -z])/π, 0.75)
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

@testset "Determine Planar Angles" begin
    @testset "2D Single Angle" begin
        @test isapprox(rot2Dangle([+1.0;  0.0]), 0.00);
        @test isapprox(rot2Dangle([+1.0; +1.0]), 0.25);
        @test isapprox(rot2Dangle([ 0.0; +1.0]), 0.50);
        @test isapprox(rot2Dangle([-1.0; +1.0]), 0.75);
        @test isapprox(rot2Dangle([-1.0;  0.0]), 1.00);
        @test isapprox(rot2Dangle([-1.0; -1.0]), 1.25);
        @test isapprox(rot2Dangle([ 0.0; -1.0]), 1.50);
        @test isapprox(rot2Dangle([+1.0; -1.0]), 1.75);
    end

    @testset "2D Angle" begin
        @test isapprox(rot2Dangle([1.;0.], [+1.; 0.]), 0.00);
        @test isapprox(rot2Dangle([1.;0.], [+1.;+1.]), 0.25);
        @test isapprox(rot2Dangle([1.;0.], [ 0.;+1.]), 0.50);
        @test isapprox(rot2Dangle([1.;0.], [-1.;+1.]), 0.75);
        @test isapprox(rot2Dangle([1.;0.], [-1.;+0.]), 1.00);
        @test isapprox(rot2Dangle([1.;0.], [-1.;-1.]), 1.25);
        @test isapprox(rot2Dangle([1.;0.], [ 0.;-1.]), 1.50);
        @test isapprox(rot2Dangle([1.;0.], [+1.;-1.]), 1.75);

        @test isapprox(rot2Dangle([1.;1.], [+1.; 1.]), 0.00);
        @test isapprox(rot2Dangle([1.;1.], [ 0.;+1.]), 0.25);
        @test isapprox(rot2Dangle([1.;1.], [-1.;+1.]), 0.50);
        @test isapprox(rot2Dangle([1.;1.], [-1.; 0.]), 0.75);
        @test isapprox(rot2Dangle([1.;1.], [-1.;-1.]), 1.00);
        @test isapprox(rot2Dangle([1.;1.], [ 0.;-1.]), 1.25);
        @test isapprox(rot2Dangle([1.;1.], [+1.;-1.]), 1.50);
        @test isapprox(rot2Dangle([1.;1.], [+1.; 0.]), 1.75);

        @test isapprox(rot2Dangle([+1.;-1.], [+1.;-1.]), 0.00);
        @test isapprox(rot2Dangle([+1.;-1.], [+1.; 0.]), 0.25);
        @test isapprox(rot2Dangle([+1.;-1.], [+1.;+1.]), 0.50);
        @test isapprox(rot2Dangle([+1.;-1.], [ 0.;+1.]), 0.75);
        @test isapprox(rot2Dangle([+1.;-1.], [-1.;+1.]), 1.00);
        @test isapprox(rot2Dangle([+1.;-1.], [-1.; 0.]), 1.25);
        @test isapprox(rot2Dangle([+1.;-1.], [-1.;-1.]), 1.50);
        @test isapprox(rot2Dangle([+1.;-1.], [ 0.;-1.]), 1.75);

        @test isapprox(rot2Dangle([-1.;+0.], [-1.; 0.]), 0.00);
        @test isapprox(rot2Dangle([-1.;+0.], [-1.;-1.]), 0.25);
        @test isapprox(rot2Dangle([-1.;+0.], [ 0.;-1.]), 0.50);
        @test isapprox(rot2Dangle([-1.;+0.], [+1.;-1.]), 0.75);
        @test isapprox(rot2Dangle([-1.;+0.], [+1.; 0.]), 1.00);
        @test isapprox(rot2Dangle([-1.;+0.], [+1.;+1.]), 1.25);
        @test isapprox(rot2Dangle([-1.;+0.], [ 0.;+1.]), 1.50);
        @test isapprox(rot2Dangle([-1.;+0.], [-1.;+1.]), 1.75);
    end
end

@testset "Point in a Square" begin
    V = [
        0.0 1.0 1.0 0.0;
        0.0 0.0 1.0 1.0
    ];
    EV = [ [1, 2], [2, 3], [3, 4], [4, 1] ];

    @test pointInPolyhedron2D([0.5;  0.5], V, EV) == +1;
    @test pointInPolyhedron2D([0.1;  0.1], V, EV) == +1;
    @test pointInPolyhedron2D([0.1;  0.9], V, EV) == +1;
    @test pointInPolyhedron2D([0.5; 1e-8], V, EV) == +1;

    @test pointInPolyhedron2D([0.5;  1.5], V, EV) == -1;
    @test pointInPolyhedron2D([1.0;  1.5], V, EV) == -1;
    @test pointInPolyhedron2D([0.5;-1e-7], V, EV) == -1;
    @test pointInPolyhedron2D([1.5;  0.5], V, EV) == -1;

    @test pointInPolyhedron2D([1.0;  0.5], V, EV) ==  0;
    @test pointInPolyhedron2D([1.0;  1.0], V, EV) ==  0;
    @test pointInPolyhedron2D([0.5; 1e-9], V, EV) ==  0;
    @test pointInPolyhedron2D([0.5;-1e-9], V, EV) ==  0;

end
