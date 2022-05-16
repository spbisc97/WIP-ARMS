function next_single_control = iLQR_function(istate, state_d, it)
    global u %use global u to remember last optimized values
    next_single_control = 0; %default next control if we have a nan value
    dt = 0.01; %delta t for integration
    t = it; %inital t
    state = istate; %inital state
    [n_states, sz] = size(state_d);

    Q = eye(n_states) * 1;
    Q(1, 1) = 20;
    Q(3, 3) = 10;
    R = 0.0001;

    iterations = 10000;
    cost = 0;
    new_cost = 0;
    j_rm = 0.001;

    horizon = 4; %time S
    horizon_disc = floor(horizon / dt);
    defects_max = ones(n_states, 1) * 0.5;
    defects_max(1, 1) = 0.2;
    defects_max(3, 1) = 0.2;
    defects_max(2, 1) = 40;
    defects_max(4, 1) = 40;

    if horizon_disc > (sz - 1)
        horizon = (sz - 1) * dt;
        horizon_disc = (sz - 1);
    end

    S = repmat(Q, 1, horizon_disc);
    L = zeros(1, horizon_disc * n_states); %size depends both from the number of controls and states
    l = zeros(1, horizon_disc); % size depends from the number of controls

    if exist("u", "var")
        [usz, ~] = size(u);
        u = [0; u];
        u = [u(3:end); zeros(horizon_disc - usz, 1)];
    else
        u = ones(horizon_disc, 1) * (0);
    end

    new_u = u;

    state_array = [istate];
    control_array = [];
    time_array = [it];

    while t < it + horizon - dt
        %compute the forward dynamics to define the defects
        elem = floor((t - it + dt) / dt);
        dy = ForwardDynamics(state, u(elem));
        state = euler_integration_fun(state, dy, dt);
        t = t + dt;
        control_array = [control_array, u(elem)];
        time_array = [time_array, t];
        state_array = [state_array, state];
    end

    %plot(time_array, state_array)
    %legend("x","dx","phi","dphi")
    %pause
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
    defects = state_array(:, 1:horizon_disc) - state_d(:, 1:horizon_disc);
    %defects(:, 1:end - 1) = 0;
    %disp("defects")
    %disp(defects(:)')
    s = zeros(n_states, horizon_disc); %deep horizon+1 and hight is n_states

    %start the optimizing iterations
    for iteration = 1:iterations - 1

        %fill the s matrix , expecially the last element
        s(:, horizon_disc) = Q * defects(:, horizon_disc); %TBC % check if it is Q*defect(horizon disc)
        %last element = horizon_disc
        A_ = zeros(n_states, (horizon_disc - 1) * n_states);
        B_ = zeros(n_states, horizon_disc - 1);
        %back
        for n = (horizon_disc - 1):-1:1
            N = (4 * n - 3):(4 * n);
            [A_(:, N), B_(:, n)] = linearization_discretization(u(n), state_array(:, n), 1); %maybe discrete
            A = A_(:, N);
            B = B_(:, n);

            %compute r,h,G,H to simplify S and s computations
            P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
            r = R * u(n); % the derivative of uRu
            q = Q * defects(:, n);
            h = r + B.' * (s(:, n + 1) ); %gu + S(:, N + 4) * defects(:, n+1)
            G = P + B.' * S(:, N + 4) * A; %Gux
            H = R + B.' * S(:, N + 4) * B; %Guu
            % disp(["H"])
            % disp([H])
            % disp("G")
            % disp([G])
            % disp("h")
            % disp([h])
            %compute Values to use in the forward iterations
            L(:, N) = -pinv(H) * G; %K
            l(:, n) = -pinv(H) * h; %d
            %disp(["L","l","n"])
            %disp([L(n),l(n),n])

            %auxiliary
            Gxx = Q + A.' * S(:, N + 4) * A;
            gx=q + A' * (s(:, n + 1) );%+ S(:, N + 4) * defects(:, n+1)

            %compute next S,s Value
            %Q+AS()A is the Gxx
            Gxu=G'
            Gux=G
            S(:, N) = Gxx + L(:, N).' * H * L(:, N) - Gxu * L(:, N) - L(:, N).' * Gux;%forse se so sbagliati con il meno

            s(:, n) =gx + G.' * l(:, n)  + L(:, N).' * (h + H * l(:, n))  ;
            % gx                                                 %gxu*d             %  %K+gu + K Guu d
            %disp(["S","s","n"])
            %disp([S(n),s(n),n])

        end

        % disp("iter")
        % disp(iteration)
        % disp("l")
        % disp(l')
        % disp("L")
        % disp(L')
        % disp("S")
        % disp(S')
        % disp("s")
        % disp(s')
        new_state_array = [istate];
        %forward iteration
        for n = 1:horizon_disc - 1
            %compute delta_u and update the value
            N = (4 * n - 3):(4 * n); %Prendere elemento da n*4-3 a elemento n*4
            A = A_(:, N);
            B = B_(:, n);
            delta_u = l(:, n) + L(:, N) * (new_state_array(:, n) - state_array(:, n)); % in this case we take the error %defects(:, n);
            new_u(n) = u(n) + delta_u;

            % new_state = state_array(:, n + 1) ... % take the previous value
            %     + (A + B * L(:, N)) * (new_state_array(:, n) - state_array(:, n)) ... %like add the control
            %     +B * l(:, n) + defects(:, n + 1);

            dy = ForwardDynamics(new_state_array(:, n), new_u(n));
            new_state = euler_integration_fun(new_state_array(:, n), dy, dt);
            new_state_array = [new_state_array, new_state];

            if (isnan(new_u(n)))
                disp(iteration)
                disp("isnan")
                %pause
                u = u * 0;
                return
                %continue
            end

        end

        t = it;
        state = istate;
        state_array = [istate];
        control_array = [];
        time_array = [it];

        while t < it + horizon -dt
            %compute the forward dynamics to define the defects
            elem = floor((t - it + dt) / dt);
            dy = ForwardDynamics(state, new_u(elem));
            state = euler_integration_fun(state, dy, dt);
            t = t + dt;
            control_array = [control_array, u(elem)];
            time_array = [time_array, t];
            state_array = [state_array, state];

        end

        %defects=distance between state array and desired state
        %(traj and desired trajectory)
        defects = state_array(:, 1:horizon_disc) - state_d(:, 1:horizon_disc);
        %defects(:, 1:end - 1) = 0;
        disp("istate")
        disp(istate)
        disp("state array")
        disp(state_array(:, :))
        disp("state_d traj")
        disp(state_d(:, 1:horizon_disc))
        disp("defects")
        disp(defects(:, :))
        disp("control")
        disp(u')
        tiledlayout(2, 1);
        nexttile
        plot(time_array, [state_array(1, :); state_array(3, :)])
        legend("x", "phi")
        nexttile
        plot(time_array, [state_array(2, :); state_array(4, :)])
        legend("dx", "dphi")

        new_cost = 0;

        for i = 1:horizon_disc - 1
            new_cost = new_cost + defects(:, i + 1)' * Q * defects(:, i + 1) + control_array(:, i)' * R * control_array(:, i);
        end

        relative = abs(new_cost - cost) / new_cost;
        disp('new_cost')
        disp(new_cost)
        disp('relative')
        disp(relative)

        pause(0.005)
        u = new_u;
        next_single_control = u(1);

        %check cost increments and return if solved
        %if (((relative < j_rm)))
        if (((relative < j_rm)) && all((abs(defects(:, horizon_disc)'))' < defects_max))

            return
        end

        cost = new_cost;

    end

    %save('ilqrVars.mat') % save variables to
end


function state = dynamics_rk4(state, u, dt)
    f1 = ForwardDynamics(state, u);
    f2 = ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end
