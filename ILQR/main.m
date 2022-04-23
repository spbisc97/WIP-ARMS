close all , clear, clc

y=1;


t=0.01;
tf=10;
dt=0.01;
y_d=t:dt:tf;
state_array=[];
control_array=[];
time_array=[];

while t<tf
    state_array = [state_array,y];
    u=LQR_function(y,y_d(floor(t/dt)));
    dy=ForwardDynamics(y,u);
    y=euler_integration_fun(y,dy,dt);
    t=t+dt;
    control_array = [control_array,u];
    time_array=[time_array,t];
end