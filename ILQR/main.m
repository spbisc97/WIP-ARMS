close all , clear, clc

y=1;


t=0.01;
tf=100;
dt=0.01;

time=t:dt:tf;
y_d=sin(time);
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
tiledlayout(1,2)
nexttile
plot(time_array,state_array)
nexttile
plot(time_array,control_array)