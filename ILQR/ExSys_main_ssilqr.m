function ExSys_main_ssilqr(Q, R, wclose)
    %if arg are less then 3 set wclose(close windows) to false
    if nargin < 4
        wclose = 0;
    end
    if nargin < 2
        R=0.1;
    end
    if nargin <1
        Q=0;
    end
    if nargin <3
        Qn=10;
        
    end

    if (wclose)
        close all
    end

    clc;
    % global u;
    u =-4;
    y = 1.5; %initial point
    t = 0.01; %initial time
    tf = 15; %final time
    dt = 0.01; %increasing time %time step

    time = t:dt:tf+dt; %time array
    %traj_d=-0.06*(time-1.7).^5+0.6; %desired trajectory
    traj_d = repmat(0, [1, (tf / dt)+1]);
    %traj_d(1,:)=-2*sin(time/2); %desired trajectory
    state_array = [y]; %array degli stati
    control_array = []; %array del controllo
    time_array = [t]; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = traj_d(:,1); %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")

    il=iLQR_GNMS(ExSys(),Q,R,Qn);
    il.order=[1];
    il.names=["x"];   
    il.plot_steps=1;  
    il.plot_start=true;
    il.plot_end=true;
    il.pieces=5;
    il.horizon=2.99;
    il.defects_max=il.defects_max*1e-6;
    il.plot_duration=0;

    while t < tf-5
        %find u control
        t_disc=floor(t / dt);
        y_des = traj_d(:, t_disc);
        u =-4;
        u = il.SS_iLQR(y,traj_d(:,t_disc:end),t,u);
        %u = LQR_function(y, y_des, Q, R);
        u_next=u(:,1);



        control_array = [control_array, u_next];

        %compute dynamics

        
        y = dynamics_rk4(y, u_next, dt);
        state_array = [state_array, y];
        y_d_array = [y_d_array, y_des];
        t = t + dt; %time increment
        time_array = [time_array, t];

    end

    save('mainVars.mat') % save variables to
    tiledlayout(3, 1)
    nexttile
    plot(time_array, state_array)
    hold on
    plot(time_array, y_d_array)
    legend("x", "dx", "phi", "dphi", "d-x", "d-dx", "d-phi", "d-dphi")
    title("trajectory and desired trajectory")

    nexttile
    plot(time_array, y_d_array - state_array)
    legend("x", "dx", "phi", "dphi")
    title("error trajectory")

    nexttile
    plot(time_array, [control_array, 0])
    title("controls")
end

function state = dynamics_rk4(state, u, dt)
    f1 = ExSys.ForwardDynamics(state, u);
    f2 = ExSys.ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ExSys.ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ExSys.ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end
