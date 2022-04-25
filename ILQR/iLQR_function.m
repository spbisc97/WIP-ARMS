function u=iLQR_function(state,state_d,it)

dt=0.01;
t=it;
y=state;

[~,sz]=size(state_d);


u=zeros(sz);


d_max=2;
J_rel=0.01;
iterations=2;
horizon=3; %time S
horizon_disc=floor(horizon/dt);
Q = 100;
R = 10;
state_array=[];
control_array=[];
time_array=[];

P_f=eye(1);
P=P_f;
if horizon_disc>sz
    horizon=floor(sz/dt);
    horizon_disc=sz;
end

while t<it+horizon
    dy=ForwardDynamics(y,u);
    y=euler_integration_fun(y,dy,dt);
    state_array = [state_array,y];
    time_array=[time_array,t];
    t=t+dt;
end
defects=state_array(1:horizon_disc)-state_d(1:horizon_disc);

  

for iteration = 1:iterations-1
    
    A_=[];
    B_=[];
    for step=1:horizon_disc
        [A_(step),B_(step)]=linearization_discretization(u(step),state_array(step));
    end
    for step=0:horizon_disc-1
        P=1;
        A=A_(horizon_disc-step);
        B=B_(horizon_disc-step);
        for i=1:3
            P_next=A'*P*A-(A'*P*B)*pinv(R+B'*P*B)*(B'*P*A)+Q; 
            P=P_next;
        end
        K=(R+B'*P*B)*(B'*P*A);
        disp('')

        u(horizon_disc-step)=-K*(state_array(horizon_disc-step)-state_d(horizon_disc-step));
    end
end
disp(u)
u=u(1);

end
