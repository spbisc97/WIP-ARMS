function u=iLQR_function(istate,state_d,it)
state_d=state_d(:);
dt=0.01;
t=it;
state=istate;
disp("ilqr")

%pause

[sz,~]=size(state_d);




Q = 10;
R = 0.01;

iterations=50;
horizon=0.06; %time S
horizon_disc=floor(horizon/dt);
if horizon_disc>sz
    horizon=sz*dt;
    horizon_disc=sz;
end

S=repmat(Q,horizon_disc,1);
L=zeros(horizon_disc,1);
l=zeros(horizon_disc,1);
u=ones(horizon_disc,1)*(0);


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
%disp("defects")
%disp(defects(:)')
s=zeros(horizon_disc,1);
s(horizon_disc)=Q*state_array(horizon_disc);

%start the optimizing iterations
for iteration = 1:iterations-1
    %pause(0.1)

    %dispvar()
    %plot(defects)
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
    for step=1:horizon_disc-1 %horizon_disc-1 times
        n=horizon_disc-step;
        P=0;%set mixed weight to zero delta_u*P*delta_x
       
        A=A_(n);
        B=B_(n);

        %compute r,h,G,H to simplify S and s computations
        r=R*u(n); %should be the derivative of uRu
        h=r+B'*S(n+1)*B;
        G=P+B'*S(n+1)*A;
        H=R+B'*(s(n+1)+S(n+1)*defects(n));
        %disp(["H","G","h"])
        %disp([H,G,h])

      

        %compute Values to use in the forward iterations
        L(n)=-pinv(H)*G;
        l(n)=-pinv(H)*h;
        %disp(["L","l","n"])
        %disp([L(n),l(n),n])

        %compute next S,s Value
        S(n)=Q+A' *S(n+1)*A-L(n)'*H*L(n);
        q=Q*defects(n);
        s(n)=q+A'*(s(n+1)+S(n+1)*defects(n))+G'*l(n)+L(n)'*(h+H*l(n));
        %disp(["S","s","n"])
        %disp([S(n),s(n),n])

    end
    disp("iter")
    disp(n)
    disp("l")
    disp(l')
    disp("L")
    disp(L')
    disp("S")
    disp(S')
    
    %forward iteration
    for n=1:horizon_disc-1
        %compute delta_u and update the value
        delta_u=l(n)+L(n)*defects(n);
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

disp("next control")
disp(u(1))
u=u(1);
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