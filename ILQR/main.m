function main(Q, R, wclose)
    %if arg are less then 3 set wclose(close windows) to false
    if nargin < 3
        wclose = 0;
    end

    if nargin < 2
        R = 0.0001;
    end

    if nargin < 1
        Q = diag([10, 1, 100, 1]);
    end

    if (wclose)
        close all
    end

    clc;
    u = 0;
    y = [0; 0; 0; 0]; %initial point
    t = 0.01; %initial time
    tf = 10; %final time
    dt = 0.01; %increasing time %time step

    time = t:dt:tf; %time array
    %traj_d=-0.06*(time-1.7).^5+0.6; %desired trajectory
    traj_d = repmat([1; 0; pi; 0], [1, tf / dt]);
    %traj_d(1,:)=-0.06*(time-1.7).^5+0.6; %desired trajectory
    state_array = [y]; %array degli stati
    control_array = []; %array del controllo
    time_array = [t]; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = [0; 0; pi; 0]; %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")

    rk4_state_array = [y];
    linearized__state_array = [y];

    while t < tf %process start

        %find u control
        y_des = traj_d(:, floor(t / dt));
        %u=iLQR_function(y,traj_d(:,floor(t/dt):end),t);
        u = iLQR_DDP_function(y, traj_d(:, floor(t / dt):end), t);
        %u = LQR_function(y, y_des, Q, R);
        %save to plot
        control_array = [control_array, u];

        %compute dynamics

        y_d_array = [y_d_array, y_des];

        dy = ForwardDynamics(y, u);
        y = euler_integration_fun(y, dy, dt);
        state_array = [state_array, y];

        %     [A,B]=linearization_discretization(u,rk4_state_array(:,end),1);
        %     linearized__state_array =[linearized__state_array, A*rk4_state_array(:,end)+B*u];
        %
        %
        %     rk4_state_array = [rk4_state_array,dynamics_rk4(rk4_state_array(:,end),u,dt)];
        %

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

    % nexttile
    % plot(time_array, rk4_state_array)
    % hold on
    % plot(time_array, y_d_array)
    % legend("x", "dx", "phi", "dphi", "d-x", "d-dx", "d-phi", "d-dphi")
    % title("trajectory and desired trajectory")

    % nexttile
    % plot(time_array, linearized__state_array)
    % hold on
    % plot(time_array, y_d_array)
    % legend("x", "dx", "phi", "dphi", "d-x", "d-dx", "d-phi", "d-dphi")
    % title("trajectory and desired trajectory")

    nexttile
    plot(time_array, y_d_array - state_array)
    legend("x", "dx", "phi", "dphi")
    title("error trajectory")

    nexttile
    plot(time_array, [control_array, 0])
    title("controls")

end

function state = dynamics_rk4(state, u, dt)
    f1 = ForwardDynamics(state, u);
    f2 = ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end
