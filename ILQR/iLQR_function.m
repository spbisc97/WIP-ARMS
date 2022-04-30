function u=iLQR_function(state,state_d,it)
state_d=state_d(:);
dt=0.01;
t=it;
y=state;
disp("ilqr")



[sz,~]=size(state_d);




Q = 100;
R = 10;
d_max=2;
J_rel=0.01;
iterations=2;
horizon=0.05; %time S
horizon_disc=floor(horizon/dt);
if horizon_disc>sz
    horizon=sz*dt;
    horizon_disc=sz;
end

S=repmat(Q,horizon_disc,1);
L=zeros(horizon_disc,1);
l=zeros(horizon_disc,1);
u=ones(horizon_disc,1)*(-0.1);

state_array=[];
control_array=[];
time_array=[];



myvals=whos;
for n = 1:length(myvals)
    if isnan(myvals(n).name)
      eval(myvals(n).name)
    end
end 

while t<it+horizon
    %compute the forward dynamics to define the defects
    elem=floor((t-it+dt)/dt);
    disp("state array")
    disp(state_array)
    disp("state")
    disp(state)
    dy=ForwardDynamics(state,u(elem));
    state=euler_integration_fun(state,dy,dt);
    control_array=[control_array;u(elem)];
    time_array=[time_array;t];
    state_array = [state_array;state];

    t=t+dt; 
end

%defects=distance between state array and desired state 
%(traj and desired trajectory)
defects=state_array(1:horizon_disc)-state_d(1:horizon_disc);

s=zeros(horizon_disc,1);
s(horizon_disc)=Q*state_array(horizon_disc);

%start the optimizing iterations
for iteration = 1:iterations-1
    dispvar()
    plot(defects)
    pause
    %fill the s matrix , expecially the last element
    %last element = horizon_disc
    

    A_=[];
    B_=[];
    %compute the linearized dynamics 
    for step=1:horizon_disc
        [A_(step),B_(step)]=linearization_discretization(u(step),state_array(step));
    end
    %backword iteration
    for step=1:horizon_disc-1 %horizon_disc-1 times
        dispvar()
        P=0;%set mixed weight to zero delta_u*P*delta_x 
        n=horizon_disc-step;
        A=A_(n);
        B=B_(n);

        %compute r,h,G,H to simplify S and s computations
        r=0; %should be the derivative of uRu
        H=R+B'*(s(n+1)+S(n+1)*defects(n));
        G=P+B'*S(n+1)*A;
        h=r+B'*S(n+1)*B;

        myvals=whos;
for ff = 1:length(myvals)
    eval(myvals(ff).name)
    if isnan(myvals(ff).name)
      eval(myvals(ff).name)
    end
end 

        %compute Values to use in the forward iterations
        L(n)=-pinv(H)*G;
        l(n)=-pinv(H)*h;
        
        %compute next S,s Value
        q=Q*defects(n);
        S(n)=Q+A'*S(n+1)*A-L(n)'*H*L(n);
        s(n)=q*A'*(s(n+1)+S(n+1)*defects(n))+G'*l(n)+L(n)'*(h+H*l(n));
    end

    %forward iteration
    for n=1:horizon_disc-1
        %compute delta_u and update the value
        delta_u=l(n)+L(n)*defects(n);
        u(n)=u(n)+delta_u;
    end
    
%     t=it;
%     y=state;
%     state_array=[];
%     control_array=[];
%     time_array=[];
% 
%     while t<it+horizon 
%     %compute the forward dynamics to define the defects
%     state_array = [state_array;y];
%     dy=ForwardDynamics(y,u);
%     y=euler_integration_fun(y,dy,dt);
%     control_array=[control_array;u];
%     time_array=[time_array,t];
%     t=t+dt;
%     end
%     %defects=distance between state array and desired state (traj
%     %and desired trajectory)
%     defects=state_array(1:horizon_disc)-state_d(1:horizon_disc);
%      


end
disp("control")
disp(u(1))
u=u(1);
dispvar()
save('ilqrVars.mat') % save variables to 

end
function dispvar()
myvals = whos;
for n = 1:length(myvals)
    if isnan(myvals(n).name)
      eval(myvals(n).name)
    end
end
end