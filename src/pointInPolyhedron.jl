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

	# Faces triangulation is performed (if needed)
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
`σ` must be a simplex, therefore, since the dimension is three, then `σ` must be
made of four vertices (non collinear).
"""
function simplexSolidAngle(v::Int64, σ::Lar.Points)::Float64
    @assert size(σ) == (3, 4)
		"ERROR: the simplex σ must be made of four 3D points.";
    p = σ[:, v];
    Σ = σ[:, setdiff(1:n, v)]
    simplexSolidAngle(p, Σ);
end

"""
    simplexSolidAngle(p::Array{Float64,1}, Σ::Lar.Points)::Float64
Evaluate the Solid Angle of `p` with respect to `Σ`.
`Σ ∪ p` must be a simplex, therefore, since 3 is the dimension, then `Σ` must be
made of `3` vertices (non collinear).

#ToDo check collinearity
"""
function simplexSolidAngle(p::Array{Float64,1}, Σ::Lar.Points)::Float64
    @assert size(Σ) == (3, 3) "ERROR: Σ must be made of three 3D points."
    @assert length(p) == 3 "ERROR: p must be a 3D point."

	# If two base vertices are equilocated then the solid angle is null.
	if |([Σ[:, i] == Σ[:, j] for i = 1:3 for j=i+1:3]...)
		return 0.;
	end

	# If `p` is equilocated with a Σ point then the solid angle is not defined.
	if |([p == Σ[:, i] for i = 1 : 3]...)
		println("The solid angle for $p with respect to $Σ is not defined.");
		return NaN;
	end

	# If `p` is coplanar with the base point then the angle is
	#  - 0 if `p` is not in the triangle defined by `Σ`
	#  - 2π otherwise
#	if areCoplanar(p, Σ[:, 1], Σ[:, 2], Σ[:, 3])
#		print("$p is coplanar with $Σ ");
#		if isInsideTriangle(p, Σ)
#			println("and it is inside of it.");
#			return 2π;
#		else
#			println("but it is outside of it.");
#			return 0.;
#		end
#	end

	# If `p` is collinear with two base point then the solid angle is null
#	if |([areCollinear(p, Σ[:, i], Σ[:, j]) for i = 1 : 3 for j = i+1 : 3]...)
#		println("$p is collinar with two points of $Σ.");
#		return 0;
#	end

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
	);
    return abs(2 * atan(tan05Ω));
end
