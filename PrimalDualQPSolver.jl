module PrimalDualQPSolver

using Convex, SCS
set_default_solver(SCSSolver(verbose=0))

function primaldual_interiorpoint_solve(Q::Array, q::Array, # assume that both Q,q are given
                                        G::Array=[], h::Array=[],
                                        A::Array=[], b::Array=[])

    # TODO add asserts, construct dummy matrices
    @assert size(Q,1) == size(Q,2)
    @assert size(Q,1) == length(q)
    # @assert length(q) == 1 # TODO figure out how to deal with n-element and nx1 arrays
    n = size(Q,1)

    # Is there a better way to do this?
    G_full = fill(0., n, n)
    h_full = fill(0., n, 1)
    if length(G)>0
        @assert size(G,2) == n
        @assert size(G,1) == length(h)
        G_full[1:size(G,1), 1:size(G,2)] = G
        h_full[1:length(h)] = h
    end

    A_full = fill(0., n, n)
    b_full = fill(0., n, 1)
    if length(A)>0
        @assert size(A,2) == n
        @assert size(A,1) == length(b)
        A_full[1:size(G,1), 1:size(G,2)] = A
        b_full[1:length(h)] = b
    end

    # Step 1: initialization
    if is_singular(A_full)
        y = zeros(n)
    else
        println("y")
        y = A_full \ b_full
    end
    if is_singular(G_full)
        z = zeros(n)
    else
        println("z")
        z = G \ (y+h_full)
    end
    if is_singular(Q)
        x = zeros(n)
    else
        println("x")
        x = Q \ (-q - G_full'*z - A_full'*y)
    end

    print(x,y,z)

    #TODO 

    # imat = eye(n)
    # zmat = zeros(n,n)
    # L_init = hvcat(3, Q, G_full', A_full', G_full, -imat, zmat, A_full, zmat, zmat)
    # println(L_init)
    # r_init = vcat(-q, h_full, b_full)
    # xzy = L_init \ r_init

    # println(svdvals(L_init))
    # println(xzy)

end

function convex_solve(Q::Array, q::Array, # assume that both Q,q are given
                 G::Array=[], h::Array=[],
                 A::Array=[], b::Array=[])
     n = size(Q,1)
     x = Variable(n)

     constraints = []
     if length(G)>0
         push!(constraints, (G*x <= h))
     end
     if length(A)>0
         push!(constraints, (A*x == b))
     end

     obj = quadform(x,Q) + q'*x
     p = minimize(obj, constraints...)
     solve!(p)

     for c in constraints
         print("dual: ")
         println(c.dual)
     end

     return (p.optval, x.value)
end

function is_singular(mat::Array)
    return minimum(svdvals(mat)) == 0
end

end #end module
