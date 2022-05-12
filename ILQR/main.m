function main(Q, R, wclose)
    %if arg are less then 3 set wclose(close windows) to false
    if nargin < 3
        wclose = 0;
    end
    if nargin < 2
        R=0.0001
    end
    if nargin <1
        Q=diag([10,1,100,1])
    end

    if (wclose)
        close all
    end

    clc;

    y = [0.5; 0; pi; 0]; %initial point
    t = 0.01; %initial time
    tf = 3; %final time
    dt = 0.01; %increasing time %time step

    time = t:dt:tf; %time array
    %traj_d=-0.06*(time-1.7).^5+0.6; %desired trajectory
    traj_d = repmat([0; 0; pi; 0], [1, tf / dt]);
    state_array = []; %array degli stati
    control_array = []; %array del controllo
    time_array = []; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = []; %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")

    while t < tf %process start

        %find u control
        y_des = traj_d(:, floor(t / dt));
        u=iLQR_function(y,traj_d(:,floor(t/dt):end),t);
        %u = LQR_function(y, y_des, Q, R);
        %save to plot
        control_array = [control_array, u];
        time_array = [time_array, t];
        y_d_array = [y_d_array, y_des];
        state_array = [state_array, y];

        %compute dynamics
        dy = ForwardDynamics(y, u);
        y = euler_integration_fun(y, dy, dt); % start to y and use a dy increment for delta t (euler integration)
        t = t + dt; %time increment
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
    plot(time_array, control_array)
    title("controls")
