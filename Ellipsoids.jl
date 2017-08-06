# Building on: https://github.com/rdeits/LoewnerJohnEllipsoids.jl/blob/master/src/LoewnerJohnEllipsoids.jl

__precompile__()
module Ellipsoids

using Convex, SCS
set_default_solver(SCSSolver(verbose=0));
using Polyhedra
using CDDLib
using PyPlot

export Ellipsoid, outer_ellipsoid, inner_ellipsoid, plot_ellipsoid, plot_points, plot_polygon

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

function outer_ellipsoid{T}(points::Vector{T}, plot::Bool=false)
    n_dims = length(points[1])
    if length(points) <= n_dims
        error("Can't make ellipsoid with <(dim+1) points")
    end

    A = Variable(n_dims, n_dims)
    b = Variable(n_dims)

    objective = logdet(A)
    problem = maximize(objective)
    for point in points
         problem.constraints += norm(A*vec(point) + b) <= 1
    end
    solve!(problem)

    A, b = A.value, b.value
    ell = Ellipsoid(A, b,  "inverse")

    if plot
        plot_points(points)
        plot_ellipsoid(ell)
    end

    return ell
end

function inner_ellipsoid{T}(points::Vector{T}, plot::Bool=false)
    n_dims = length(points[1])
    if length(points) <= n_dims
        error("Can't make ellipsoid with <(dim+1) points")
    end

    # Convert points to inequality form
    x = list_to_stack(points)
    v_rep = SimpleVRepresentation(x)
    poly = polyhedron(v_rep, CDDLibrary())
    h_rep = SimpleHRepresentation(hrep(poly))
    Ah = h_rep.A
    bh = h_rep.b

    if plot
        plot_polygon(points)
    end

    return inner_ellipsoid(Ah, bh, plot)
end

function inner_ellipsoid{T}(A::Array{T,2}, b::Array{T,1}, plot::Bool=false)
    n_dims = size(A)[2]
    B = Variable(n_dims,n_dims)
    d = Variable(n_dims)

    objective = logdet(B)
    problem = maximize(objective)
    for i in 1:size(A)[1]
         problem.constraints += norm(B*A[i,:])+A[i,:]'*d <= b[i]
    end
    solve!(problem)

    B,d = B.value,d.value
    ell = Ellipsoid(B, d, "forward")

    if plot
        plot_ellipsoid(ell)
    end

    return ell
end

# Use forward image representation to plot, since it's the easiest to understand
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
                if isempty(u_vecs)
                    u_vecs = v
                end
                u_vecs = hcat(u_vecs, v)
            end
        end
    else
        println("Can't plot ellipsoids that aren't 2D or 3D")
        return
    end
     e = A*u_vecs + repmat(c, 1, size(u_vecs,2))

     if n_dims==2
         plot(e[1,:], e[2,:], "b.")
     elseif n_dims==3
         plot3D(e[1,:], e[2,:], e[3,:], "b.")
     end
     axis("equal")
end

function plot_points{T}(points::Vector{T})
    n_dims = length(points[1])
    x = list_to_stack(points)
    if n_dims==2
        plot(x[:,1], x[:,2], "r.")
    elseif n_dims==3
        plot3D(x[:,1], x[:,2], x[:,3], "r.")
    end
end

function plot_polygon{T}(points::Vector{T})
    n_dims = length(points[1])

    x = list_to_stack(points)

    if n_dims==2
        xplot = [x[:,1]; x[1,1]; x[end,1]]
        yplot = [x[:,2]; x[1,2]; x[end,2]]
        plot(xplot, yplot, "r")
    elseif n_dims==3
        xplot = [x[:,1]; x[1,1]; x[end,1]]
        yplot = [x[:,2]; x[1,2]; x[end,2]]
        zplot = [x[:,3]; x[1,3]; x[end,3]]
        plot3D(xplot, yplot, zplot, "r.")
    end
end

# TODO do this better
function list_to_stack{T}(v::Vector{T})
    n_dims = length(v[1])
    stack = zeros(length(v), n_dims)
    for i in 1:length(v)
        for j in 1:n_dims
            stack[i, j] = v[i][j]
        end
    end
    return stack
end

# TODO run instantiation of Float64 to precompile?

end # module
