

V, (VV,EV,FV,CV) = Lar.simplex(3, true);
#==
tetra = V, EV,FV,CV;
twotetra = Lar.Struct([ tetra, Lar.t(1.,1.,1.), tetra ]);
V,EV,FV,CV = Lar.struct2lar(twotetra);
==#
#==
cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);

triangulated_faces = Lar.triangulate(convert(Lar.Points, V'), [cop_EV, cop_FE])

FE = Array{Int64,1}[
    [1, 2, 3],
    [1, 4, 5],
    [2, 4, 6],
    [3, 5, 6, 7, 8, 9],
    [7, 8, 9],
    [7, 10, 12],
    [7, 11, 13, 14],
    [8, 10, 15],
    [8, 11, 16, 17],
    [9, 12, 15],
    [9, 13, 16, 18],
    [14, 17, 18]
];

triFV = Array{Array{Int64,1},1}[
    [[2, 3, 1]], [[4, 2, 1]], [[4, 3, 1]],
    [[2, 3, 5]], [[6, 3, 7]], [[5, 3, 6]],
    [[7, 4, 2]], [[4, 7, 3]], [[7, 2, 5]],
    [[5, 9, 10]], [[10, 6, 5]], [[11, 7, 9]],
    [[5, 9, 7]], [[11, 7, 10]], [[6, 10, 7]],
    [[9, 10, 11]]
];
myV = [
    0.0 1.0 0.0 0.0 0.5 0.25 0.25 0.25 1.25 0.25 0.25;
    0.0 0.0 1.0 0.0 0.25 0.5 0.25 0.25 0.25 1.25 0.25;
    0.0 0.0 0.0 1.0 0.25 0.25 0.5 0.25 0.25 0.25 1.25
];
GL.VIEW(GL.GLExplode(myV,myFV,1.5,1.5,1.5,99));
==#

p = [.1;.1;.1]
