function u_next=iLQR_function(istate,state_d,it)
global u
state_d=state_d(:);
dt=0.01;
t=it;
state=istate;
%disp("ilqr")

%pause

[sz,~]=size(state_d);




Q = 10;
R = 0.01;


iterations=20;
horizon=0.1; %time S
horizon_disc=floor(horizon/dt);
if horizon_disc>sz
    horizon=sz*dt;
    horizon_disc=sz;
end

S=repmat(Q,horizon_disc+1,1);
L=zeros(horizon_disc,1);
l=zeros(horizon_disc,1);
if exist("u","var")
    [usz,~]=size(u);
    u=[0;u];
    u=[u(3:end);zeros(horizon_disc-usz+1,1)];
else
    u=ones(horizon_disc,1)*(0);
end

state_array=[];
control_array=[];
time_array=[];



% myvals=whos;
% for n = 1:length(myvals)
%     if isnan(myvals(n).name)
%       eval(myvals(n).name)
%     end
% end

while t<it+horizon
    %compute the forward dynamics to define the defects
    elem=floor((t-it+dt)/dt);
    dy=ForwardDynamics(state,u(elem));
    state=euler_integration_fun(state,dy,dt);
    control_array=[control_array;u(elem)];
    time_array=[time_array;t];
    state_array = [state_array;state];
    t=t+dt;
end
% disp("istate")
% disp(state')
% disp("state array")
% disp(state_array')
% disp("state_d traj")
% disp(state_d(1:horizon_disc)')
% 
% plot(time_array,state_array)


%defects=distance between state array and desired state
%(traj and desired trajectory)
defects=state_array(1:horizon_disc)-state_d(1:horizon_disc);
% defects(1:end-1)=0;
%disp("defects")
%disp(defects(:)')
s=zeros(horizon_disc+1,1);
s(horizon_disc+1,1)=Q*state_array(horizon_disc,1);

%start the optimizing iterations
for iteration = 1:iterations-1
    %pause(0.1)

    %dispvar()
    %plot(defects)cos(y(3))
    %pause
    %fill the s matrix , expecially the last element
    %last element = horizon_disc
    A_=[];
    B_=[];
    %compute the linearized dynamics
    for step=1:horizon_disc
        [A_(step),B_(step)]=linearization_discretization(u(step),state_array(step));
    end
    %backword iteration
    for step=1:horizon_disc %horizon_disc-1 times
        n=horizon_disc-step+1;
        P=0;%set mixed weight to zero delta_u*P*delta_x
       
        A=A_(n);
        B=B_(n);

        %compute r,h,G,H to simplify S and s computations
        r=R*u(n); %should be the derivative of uRu
        h=r+B'*(s(n+1)+S(n+1)*defects(:,n));
        G=P+B'*S(n+1)*A;
        H=R+B'*S(n+1)*B;
        %disp(["H","G","h"])
        %disp([H,G,h])

      

        %compute Values to use in the forward iterations
        L(n)=-pinv(H)*G;
        l(n)=-pinv(H)*h;
        %disp(["L","l","n"])
        %disp([L(n),l(n),n])

        %compute next S,s Value
        S(n)=Q+A'*S(n+1)*A-L(n)'*H*L(n);
        q=Q*defects(:,n);
        s(n)=q+A'*(s(n+1)+S(n+1)*defects(:,n))+G'*l(n)+L(n)'*(h+H*l(n));
        %disp(["S","s","n"])
        %disp([S(n),s(n),n])

    end
    disp("iter")
    disp(iteration)
    disp("l")
    disp(l')
    disp("L")
    disp(L')
    disp("S")
    disp(S')
    disp("s")
    disp(s')
    
    %forward iteration
    for n=1:horizon_disc
        %compute delta_u and update the value
        delta_u=l(n)+L(n)*defects(:,n));
        u(n)=u(n)+delta_u;
    end

    t=it;
    state=istate;
    state_array=[];
    control_array=[];
    time_array=[];


    while t<it+horizon
        %compute the forward dynamics to define the defects
        elem=floor((t-it+dt)/dt);
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
%     defects(1:end-1)=0;
    disp("istate")
    disp(istate')
    disp("state array")
    disp(state_array')
    disp("state_d traj")
    disp(state_d(1:horizon_disc)')
    disp("defects")
    disp(defects(:)')
    disp("control")
    disp(u')
    plot(time_array,state_array)
    pause(0.1)



end

% disp("next control")
% disp(u(1))
u_next=u(1);
%save('ilqrVars.mat') % save variables to

end



function dispvar()
myvals = whos;
for n = 1:length(myvals)
    if isnan(myvals(n).name)
        eval(myvals(n).name)
    end
end
end