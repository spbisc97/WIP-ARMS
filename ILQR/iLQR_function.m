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
S=repmat(Q,horizon_disc,1);
L=repmat(0,horizon_disc,1);
l=repmat(0,horizon_disc,1);

state_array=[];
control_array=[];
time_array=[];horizon

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
    s=repmat(S.*state_array);


    A_=[];
    B_=[];
    for step=1:horizon_disc
        [A_(step),B_(step)]=linearization_discretization(u(step),state_array(step));
    end
    for step=1:horizon_disc
        P=0;%set mixed weight to zero
        n=horizon_disc-step;
        A=A_(n);
        B=B_(n);
        r=R*state_array(n);
        H=R+B'*(s(n+1)+S(n+1)*defects(n));
        G=P+B'*S(n+1)*A;
        h=r+B'*S(n+1)*B;

        L(n)=-pinv(H)*h;
        l(n)=-pinv(H)*G;

        q=Q*state_array(n);
        S(n)=Q+A'*S(n+1)*A-L'*H*L;
        s(n)=q*A'*(s(n+1)+S(n+1)*defects(n))+G'*l(n)+L(n)'*(h+H*l);
    end
    for step=1:horizon_disc
        
    
    end


end
disp(u)
u=u(1);

end
