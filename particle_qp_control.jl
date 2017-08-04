using JuMP, Ipopt
using Plots
gr()

# Constants
Δt = 0.1
num_time_steps = 20
max_acceleration = 0.5
kp = 1.5
kd = 2*sqrt(1.)

function qp_controller(pos_cur, vel_cur, pos_des, vel_des)
    model = Model(solver=IpoptSolver(print_level=0))

    # Decision variables
    @variables model begin
        -max_acceleration <= acceleration[1:2] <= max_acceleration
    end

    # Manual warm-start
#     if isdefined(:prev_solution)
#         setvalue(acceleration, prev_solution)
#     end

    # Don't really need dynamics constraint because the dynamics model is just F=m*a

    # Cost function: spring damper behavior (PD controller) on \ddot{p} between current and goal
    # TODO try converting this to Q, C matrices
    @objective(model, Min, sum((kp*(pos_des - pos_cur) + kd*(vel_des - vel_cur) - acceleration).^2))

    solve(model)
#     global prev_solution = getvalue(acceleration)
    return getvalue(acceleration)
end

pos_des = [-1,1]
vel_des = [0,0]
q = [1,0]
v = [0,0]

anim = @animate for i in 1:100
    # Plot the current and goal positions
    plot([pos_des[1]], [pos_des[2]], marker=(:circle, 10), xlim=(-2.1, 2.1), ylim=(-2.1, 2.1), legend=:none)
    plot!([q[1]], [q[2]], marker=(:hex, 10))

    # Run the MPC control optimization
#     tic()
    accel = qp_controller(q, v, pos_des, vel_des)
#     print(i, " ")
#     toc()

    # Apply the planned acceleration and simulate one step in time
    v += accel * Δt
    q += v * Δt
end

gif(anim, "img/qp_traj.gif")
