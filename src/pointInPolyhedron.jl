## Point in Polyhedron Algorithm
#  By Lane, Magedson and Rarick

using SparseArrays;
using LinearAlgebraicRepresentation, ViewerGL;
Lar = LinearAlgebraicRepresentation;
GL =  ViewerGL;

"""
	pointInPolyhedron3D(
		p::Array{Float64,1},
		V::Lar.Points,
		EV::Lar.Cells,
		FV::Lar.Cells
	)::Int64
Evaluate if `p` belongs to the complex defined by `[V, EV, FV]`.

This method computes if a given point `p` is inner, outer or edger with respect
to the complex defined by the sets `V, EV, FV`. Respectivelly it return
the integer value `+1`, `-1` or `0`.

_Note._ If something weired happens then `-3` is returned and the value for the
solid angle related to `p` is shown.
"""
function pointInPolyhedron3D(
		p::Array{Float64,1},
		V::Lar.Points,
		EV::Lar.Cells,
		FV::Lar.Cells
	)::Int64

	# Faces triangulation (if needed)
	if max([length(σ) for σ in FV]...) > 3
		cop_EV = Lar.coboundary_0(EV::Lar.Cells);
		cop_EW = convert(Lar.ChainOp, cop_EV);
		cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    	triFV = Lar.triangulate(convert(Lar.Points, V'), [cop_EV, cop_FE]);
	else
		triFV = [[σ] for σ in FV]
	end

	# Solid Angle Initialization
	Ω = 0.0;

	# For the triangulation of each face ...
    for face_complex in triFV
		# ... compute the axis face
        axis = Lar.cross(
			V[:, face_complex[1][2]] .- V[:, face_complex[1][1]],
			V[:, face_complex[1][3]] .- V[:, face_complex[1][1]]
		);
		# If the point lies on a face plan ...
		if Lar.dot(axis, V[:, face_complex[1][1]]) == Lar.dot(axis, p)
			# ... If for some face ...
			face_ref = ∪(face_complex...);
			for Σ in face_complex
				position = pointInPolyhedron2D(V[:, Σ])
				# ... the point is inner with respect to that face
				#     than it is inner to the polyhedron;
				if position == 1
					return (0, face_ref);
				# ... if it is on an edge of that polyhedron
				#     than it is either inner to that face or to that edge;
				elseif position ∋ 0;
					if position[2] ∈ EV
						return position;
					else
						return (0, face_ref);
					end
				end
			end
			# ... however, if no faces contains it, then the contribute of
			#     that face to the solid angle is equal to zero.
		# ... else if the point is not coplanar to a face ...
		else
			# ... the normal to the face is used to determine if `p` is inner
			#     (i.e. it is opposite to the side of the axis) or it is outer
			isInner = Lar.dot(axis, V[:, face_complex[1][1]]) < Lar.dot(axis, p);
        	for Σ in face_complex
				@show isInner (-1)^!isInner
				@show simplexSolidAngle(p, V[:, Σ])/π
	            Ω = Ω + (-1)^isInner * simplexSolidAngle(p, V[:, Σ])/π;
	        end
		end
    end
    if isapprox(abs(Ω), 4.0; atol=1e-1)
        return 1;
    elseif isapprox(abs(Ω), 2.0; atol=1e-1)
        return 0;
    elseif isapprox(abs(Ω), 0.0; atol=1e-1)
        return -1;
    else
        @show Ω;
		return -3;
    end

end

"""
    simplexSolidAngle(v::Int64, σ::Lar.Points)::Float64
Evaluate the Solid Angle for the `v`-th vertex of the simplex `σ`.
`σ` must be a simplex, therefore if `dim` is the dimension, then `σ` must be
made of `dim+1` vertices (non collinear).
"""
function simplexSolidAngle(v::Int64, σ::Lar.Points)::Float64
    dim, n = size(σ);
    @assert n == dim + 1 "ERROR: σ do not define a simplex."
    p = σ[:, v];
    Σ = σ[:, setdiff(1:n, v)]
    simplexSolidAngle(p, Σ);
end

"""
    simplexSolidAngle(p::Array{Float64,1}, Σ::Lar.Points)::Float64
Evaluate the Solid Angle of `p` with respect to `Σ`.
`Σ ∪ p` must be a simplex, therefore if `dim` is the dimension, then `Σ` must be
made of `dim` vertices (non collinear).

#ToDo check collinearity
#ToDo implement Lar.Cross for dimensions different from 3.
"""
function simplexSolidAngle(p::Array{Float64,1}, Σ::Lar.Points)::Float64
    dim, n = size(Σ);
    @assert n == dim "ERROR: σ do not define a simplex."
    @assert length(p) == dim "ERROR: p has a different dimension than Σ."
	# Evaluate vectors from p to Σ vertices
    A = Σ[:, 1] .- p;
    B = Σ[:, 2] .- p;
    C = Σ[:, 3] .- p;
	# Evaluate magnitudes
	a = Lar.norm(A);
	b = Lar.norm(B);
	c = Lar.norm(C);
	tan05Ω = Lar.dot(A, Lar.cross(B, C))/(
		a*b*c + Lar.dot(A, B)*c + Lar.dot(B, C)*a + Lar.dot(C, A)*b
	)
    return abs(2 * atan(tan05Ω));
end

#==
function simplexSolidAngle(p::Array{Float64,1}, Σ::Lar.Points)::Float64
    dim, n = size(Σ);
    @assert n == dim "ERROR: σ do not define a simplex."
    @assert length(p) == dim "ERROR: p has a different dimension than Σ."
    # Evaluate vectors from p to Σ vertices
    #  pa = Σ[:, 1] .- p;
    #  pb = Σ[:, 2] .- p;
    #  pc = Σ[:, 3] .- p;
    Σvec = Σ .- p;
    # Evaluate faces normal
    #  ηab = Lar.cross(pa, pb);
    #  ηbc = Lar.cross(pb, pc);
    #  ηca = Lar.cross(pc, pa);
    # do note that i%3+1 maps 1 ↦ 2, 2 ↦ 3, ..., n-1 ↦ n, n ↦ 1.
    η = [Lar.cross(Σvec[:, i], Σvec[:, i%dim+1]) for i = 1 : dim];
    # Evaluate angles between normals
    #  θa = acos(Lar.dot(ηca, ηab)/(Lar.norm(ηca)*Lar.norm(ηab)));
    #  θb = acos(Lar.dot(ηab, ηbc)/(Lar.norm(ηab)*Lar.norm(ηbc)));
    #  θc = acos(Lar.dot(ηbc, ηca)/(Lar.norm(ηbc)*Lar.norm(ηca)));
    # do note that (i+1)%n+1 maps 1 ↦ n, 2 ↦ 1, ..., n ↦ n-1.
	cosθ = [
        Lar.dot(η[(i+1)%dim+1], η[i])/
        (Lar.norm(η[(i+1)%dim+1])*Lar.norm(η[i]))
        for i = 1 : dim
    ];
	θ = acos.(cosθ);
	@show η cosθ θ
    # Return the solid angle
    #  return θa + θb + θc - π;
    return sum(θ) - π;
end
==#
