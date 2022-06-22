clc
clear all

ts=0.001;
y_i = [0, 0, 0, 0, 0, 0];
tspan = 0:ts:10;
u_i = [0; 0];

disp('linearization in 0')

discrete = false;
[A, B] = Twip.linearization_discretization(u_i, y_i, discrete);

Q = diag([10, 1, 10, 1, 10, 1]);
R = diag([1, 1]);

k = lqr(A, B, Q, R);
start = [0.1; 0; 0; 0; 0; 0];

[time, ode_y] = ode45(@( t,y)Twip.ForwardDynamics( y, -k * y), tspan, start);
y=start;
stat=[y];
for t=1:1:length(tspan)-1
    dy=Twip.ForwardDynamics( y, -k * y);
    y=Twip.euler_integration_fun(y,dy,ts);
    stat=[stat,y];

end
y=start;
sta=[y];

for t=1:1:length(tspan)-1
    
    [A,B]=Twip.linearization_discretization( -k * y,y,1);
    y_lin=A*y+B*(-k * y);
    sta=[sta,y_lin];
    if mod(t,1)~=0
        y=y_lin
    else
    dy=Twip.ForwardDynamics( y, -k * y);
    y=Twip.euler_integration_fun(y,dy,ts);
    end
end
y=start;
sta_l=[y];

% for t=1:1:1000
%     [A,B]=Twip.linearization_discretization( -k*y,y,0);
%     k = lqr(A, B, Q, R);
%     dy=Twip.ForwardDynamics( y, -k * y);
%     y=Twip.euler_integration_fun(y,dy,ts);
%     sta_l=[sta_l,y];    
% end
tiledlayout(4,1)
nexttile
plot(time, ode_y )
length(ode_y)
legend("phi", "phi_dot", "x", "x_dot", "theta", "theta_dot")
nexttile
length(stat)

plot(time, stat)
legend("phi", "phi_dot", "x", "x_dot", "theta", "theta_dot")

nexttile
length(sta)

plot(time, sta)
legend("phi", "phi_dot", "x", "x_dot", "theta", "theta_dot")
% nexttile
% length(sta_l)
% 
% plot(time, sta_l)
% legend("phi", "phi_dot", "x", "x_dot", "theta", "theta_dot")