using PyPlot

A = [3 0; 0 1]
b = [0 0]

angles = collect(0 : 0.1 : 2*pi)
vecs = [cos.(angles)-b[1] sin.(angles)-b[2]]'
ellipse = inv(A) * vecs

plot(ellipse[1,:], ellipse[2,:], "b.")
axis("equal")
