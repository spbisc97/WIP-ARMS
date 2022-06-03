using LinearAlgebra
using PyPlot
using ForwardDiff
using RobotZoo
using RobotDynamics


#AcrobotDynamics
a = RobotZoo.Acrobot()
h = 0.05


function dynamics_rk4(x, u)
    #rk4 integration with zero-order hold on u
    f1 = dynamics(a, x, u)
    f2 = dynamics(a, x + 0.5 * h * f1, u)
    f3 = dynamics(a, x + 0.5 * h * f2, u)
    f4 = dynamics(a, x + h * f3, u)
    return x + (h / 6.0) * (f1 + 2 * f2 + 2 * f3 + f4)
end

Nx = 4 #number of states
Nu = 1 #number of control 
Tfinal = 5.0
Nt = Int(Tfinal / h) + 1
thist = Array(range(0, h * (Nt - 1), step=h));



#Cost wheights 
Q = Diagonal([1.0 * ones(2); 0.1 * ones(2)]);
R = 0.001;
Qn = Array(10000.0 * I(Nx));



function stage_cost(x, u)
    return 0.5 * ((x - xgoal)' * Q * (x - xgoal)) + 0.5 * R * u * u
end

function terminal_cost(x)
    return 0.5 * ((x - xgoal)' * Qn * (x - xgoal))
end

function cost(xtraj, utraj)
    J = 0.0
    for k = 1:(Nt-1)
        J += stage_cost(xtraj[:, k], utraj[k])
    end
    J += terminal_cost(xtraj[:, Nt])
    return J
end


#initial guess

x0 = [-pi / 2; 0; 0; 0]
xgoal = [pi / 2; 0; 0; 0]
xtraj = kron(ones(1, Nt), x0)
utraj = randn(Nt - 1)

#initial rollout
for k = 1:(Nt-1)
    xtraj[:, k+1] .= dynamics_rk4(xtraj[:, k], utraj[k])
end

J = cost(xtraj, utraj)

#DDP Algo

p = zeros(Nx, Nt)
P = zeros(Nx, Nx, Nt)
d = ones(Nt - 1)
K = zeros(Nu, Nx, Nt - 1)
DJ = 0.0

gx = zeros(Nx)
gu = 0.0
Gxx = zeros(Nx, Nx)
Guu = 0.0
Gxu = zeros(Nx)
Gux = zeros(Nx)



iter = 0
while (maximum(abs.(d[:])) > 1e-3)
    global iter, J, p, P, d, K, DJ, gx, gu, Gxx, Gux, Gxu, Guu,Qn
    
    iter += 1

    p = zeros(Nx, Nt)
    P = zeros(Nx, Nx, Nt)
    d = ones(Nt - 1)
    K = zeros(Nu, Nx, Nt - 1)
    DJ = 0.0

    p[:, Nt] = Qn * (xtraj[:, Nt] - xgoal)
    P[:, :, Nt] = Qn


    #Backward pass
    for k = (Nt-1):-1:1
        #Calculate derivatives 
        q = Q * (xtraj[:, k] - xgoal)
        r = R * utraj[k]

        A = ForwardDiff.jacobian(dx -> dynamics_rk4(dx, utraj[k]), xtraj[:, k])
        B = ForwardDiff.derivative(du -> dynamics_rk4(xtraj[:, k], du), utraj[k])

        gx = q + A' * p[:, k+1]
        gu = r + B' * p[:, k+1]

        Gxx = Q + A' * P[:, :, k+1] * A
        Guu = R + B' * P[:, :, k+1] * B
        Gxu = A' * P[:, :, k+1] * B
        Gux = B' * P[:, :, k+1] * A

        d[k] = Guu \ gu
        K[:, :, k] .= Guu \ Gux

        p[:, k] .= dropdims(gx - (K[:, :, k]' * gu) + K[:, :, k]' * Guu * d[k] - Gxu * d[k], dims=2)

        P[:, :, k] .= Gxx + K[:, :, k]' * Guu * K[:, :, k] - Gxu * K[:, :, k] - K[:, :, k]' * Gux

        DJ += gu' * d[k]
    end

    #Forward
    if any(isnan,K[:, :, 1])
        global utraj
        utraj = randn(Nt - 1)
        for k = 1:(Nt-1)
            xtraj[:, k+1] .= dynamics_rk4(xtraj[:, k], utraj[k])
        end
        d[:].=1
        continue
    end
    xn = zeros(Nx, Nt)
    xn[:, 1] = xtraj[:, 1]
    un = zeros(Nt - 1)
    α = 0.5
    found_sol = false
    while !found_sol
        for k = 1:1:(Nt-1)
            un[k]= utraj[k] - α * d[k] - dot(K[:, :, k], xn[:, k] - xtraj[:, k])
            try
                xn[:, k+1] .= dynamics_rk4(xn[:, k], un[k])
            catch Exception
                un[Nt-1] = NaN
                display("catch")
                sleep(0.001)
                break
            end
            

        end
        display("α")
        display(α)
        display(K[:, :, 1])
        
        if isnan(un[Nt-1])
            α = α / 2
        else
            found_sol = true
        end
    end

    Jn = cost(xn, un)


    # while Jn > (J - 1e-2 * α * DJ)
    #     α = 0.5 * α
    #     for k = 1:(Nt-1)
    #         un[k] = utraj[k] - d[k] + α * dot(K[:, :, k], xn[:, k] - xtraj[:, k])
    #         xn[:, k+1] = dynamics_rk4(xn[:, k], un[k])
    #     end
    #     Jn = cost(xn, un)
    # end
    display(maximum(abs.(d[:])))

    J = Jn
    xtraj .= xn
    utraj .= un
end

using TrajOptPlots, Colors
using MeshCat
using StaticArrays
vis = Visualizer()
render(vis)


TrajOptPlots.set_mesh!(vis, a)
X1 = [SVector{4}(x) for x in eachcol(xtraj)]
visualize!(vis, a, thist[end], X1)
render(vis)


