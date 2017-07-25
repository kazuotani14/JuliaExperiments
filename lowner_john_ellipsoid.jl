# Finds minimum volume ellipsoid around a finite set of points (2d in this case)

# An ellipse can be represented as {v | ||Av+b||_2 <= 1}
# Constraint to contain all points is: ||A*x_i + b ||_2 <= 1 for all x_i

# We can minimize the volume of the ellipse -> minimize logdet(inv(A))
# Volume of ellipse is proportional to det(inv(A))
#    Intuition for this is that as values in A get larger, v must be smaller to keep norm <=1
# But we can't do inv in Convex so we instead maximize logdet(A)
# logdet also adds an implicit constraint that A is symmetric positive-definite

using PyPlot
using Convex, SCS
set_default_solver(SCSSolver(verbose=0))

n = 3

# Points
if n==2
     x = [10 0;
          0 5;
          1 1]
elseif n==3
     x = [10 0 0;
          0 10 0
          0 0 10;
          10 10 10]
else
     x = rand(n+1, n)
end

# Set up problem
A = Variable(n,n)
b = Variable(n)

objective = logdet(A)

constraints = (norm(A*x[1,:]+b) <= 1)
for i in 2:size(x)[1]
     constraints += (norm(A*x[i,:]+b) <= 1)
end

problem = maximize(objective, constraints)

solve!(problem)

# Results
println("Problem status: ", problem.status)
println("Optimal value: ", problem.optval)
println("A: ", A.value)
println("b: ", b.value)
A = A.value
b = b.value

# Plot

angles = collect(0 : 0.1 : 2*pi)
if n==2
     vecs = [cos.(angles)-b[1] sin.(angles)-b[2]]'
     ellipse = inv(A) * vecs

     for i in 1:size(x)[1]
          plot(x[i,1], x[i,2], "ro")
     end

     plot(ellipse[1,:], ellipse[2,:], "b.")
     axis("equal")
elseif n==3
     vecs = []
     for i in angles
          for j in angles
               v = [cos.(i).*sin.(j)-b[1] sin.(i).*sin.(j)-b[2] cos.(j)-b[3]]
               if isempty(vecs)
                    vecs = v
               end
               vecs = [vecs; v]
          end
     end

     for i in 1:size(x)[1]
          plot3D([x[i,1]], [x[i,2]], [x[i,3]], "ro")
     end

     ellipse = inv(A)*vecs'
     plot3D(ellipse[1,:], ellipse[2,:], ellipse[3,:], "b.", markersize=1)
     axis("equal")
end
