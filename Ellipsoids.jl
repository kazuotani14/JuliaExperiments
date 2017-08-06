# Building on: https://github.com/rdeits/LoewnerJohnEllipsoids.jl/blob/master/src/LoewnerJohnEllipsoids.jl

# __precompile__()
module Ell

    using Convex, SCS
    set_default_solver(SCSSolver(verbose=0));
    using Polyhedra
    using CDDLib
    using PyPlot

    export Ellipsoid

    type Ellipsoid{T}
        # Represents an ellipsoid with L (cholesky factorization of quadratic form matrix), and c
        # See https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf for details
        L::Array{T, 2}
        c::Array{T, 1}
        n_dims::Int64
    end

    Ellipsoid{T}(A::Array{T, 2}, b::Array{T, 1}, rep::String) = begin
        @assert size(A,2) == length(b)

        if rep=="forward"
            L = inv(A)'
            c = b
        elseif rep=="inverse"
            L = A'
            c = -inv(L')*b
        elseif rep=="quadratic"
            L = Vector(chol(A))'
            c = b
        else
            error("Representation must be one of: forward, inverse, quadratic")
        end

        n_dims = size(A,2)
        Ellipsoid(L, c, n_dims)
    end

    Ellipsoid{T}(A::Array{T, 2}, c::Array{T, 2}, rep::String) = Ellipsoid(A, vec(c), rep)

    # # Use forward image representation to plot, since it's the easiest to understand
    function plot_ellipsoid{T}(E::Ellipsoid{T})
        A = inv(E.L)'
        c = E.c

        n_dims = size(A,2)

        if n_dims==2
            angles = collect(0 : 0.3 : 2*pi)
            u_vecs = [cos.(angles) sin.(angles)]'
        elseif n_dims==3
            angles = collect(0 : 0.3 : 2*pi)

            # TODO make this concise...
            u_vecs = []
            for i in angles
                for j in angles
                    v = [cos.(i).*sin.(j); sin.(i).*sin.(j); cos.(j)]
                    if isempty(vecs)
                        u_vecs = v
                    end
                    u_vecs = hcat(vecs, v)
                end
            end
        else
            println("Can't plot ellipsoids that aren't 2D or 3D")
            return
        end

         e = A*u_vecs + repmat(c, 1, length(angles))

         if n_dims==2
             plot(e[1,:], e[2,:], "b.")
         elseif n_dims==3
             plot3D(e[1,:], e[2,:], e[3,:], "b.")
         end
         axis("equal")
    end

    function outer_ellipsoid{T}(points::Vector{T})
        @assert length(points) > 1

        n_dims = length(points[1])

        A = Variable(n_dims, n_dims)
        b = Variable(n_dims)

        objective = logdet(A)
        problem = maximize(objective)

        for point in points
             problem.constraints += norm(A*point + b) <= 1
        end
        solve!(problem)

        A, b = A.value, b.value

        c = A \ b # TODO check this

        return Ellipsoid(A, c,  "inverse")
    end

    # TODO run instantiation of Float64 to precompile

end # module
