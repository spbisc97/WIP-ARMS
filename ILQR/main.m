function main(Q, R, wclose)
    %if arg are less then 3 set wclose(close windows) to false
    if nargin < 3
        wclose = 0;
    end
    if nargin < 2
        R=0.0001;
    end
    if nargin <1
        Q=diag([10,1,100,1]);
    end

    if (wclose)
        close all
    end

    clc;
    % global u;
    u = [0,0];
    y = [0; 0; pi+0.1; 0]; %initial point
    t = 0.01; %initial time
    tf = 15; %final time
    dt = 0.01; %increasing time %time step

    time = t:dt:tf+dt; %time array
    %traj_d=-0.06*(time-1.7).^5+0.6; %desired trajectory
    traj_d = repmat([1; 0; pi; 0], [1, (tf / dt)+1]);
    %traj_d(1,:)=-sin(time/2); %desired trajectory
    state_array = [y]; %array degli stati
    control_array = []; %array del controllo
    time_array = [t]; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = traj_d(:,1); %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")
    while t < tf-5
        %find u control
        t_disc=floor(t / dt);
        y_des = traj_d(:, t_disc);
        u=iLQR_function(y,traj_d(:,t_disc:end),t,u);
        %u = iLQR_DDP_function(y, traj_d(:, floor(t / dt):end), t);
        %u = LQR_function(y, y_des, Q, R);
        %save to plot
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


    y = [0; 0; pi; 0]; %initial point
    t = 0.01; %initial time
    state_array = [y]; %array degli stati
    control_array = []; %array del controllo
    time_array = [t]; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = traj_d(:,1); %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")
end

function state = dynamics_rk4(state, u, dt)
    f1 = ForwardDynamics(state, u);
    f2 = ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end
