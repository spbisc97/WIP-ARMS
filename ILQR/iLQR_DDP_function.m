function u = iLQR_DDP_function(istate, state_d, it,u)
    %global  u new_u Gxx Gxu Guu Gux gx gu d L A_ A B_ B N n
    dt = 0.01; %delta t for integration

    [n_states, sz] = size(state_d);
    n_controls = 1;
    Q = eye(n_states) * 0.01;
    Q(1, 1) = 1;
    Q(3, 3) = 1;
    R = 0.001;
    alpha=0.5;

    horizon = 10; %time S
    horizon_disc = floor(horizon / dt);

    if horizon_disc > (sz) %controllo se lo stato desiderato 
        % horizon = (sz - 1) * dt;
        horizon_disc = (sz);
    end

    state_d = state_d(:, 1:horizon_disc);

    S = zeros(n_states, horizon_disc * n_states);
    s = zeros(n_states, horizon_disc);

    L = ones(n_controls, (horizon_disc - 1) * n_states); %size depends both from the number of controls and states
    d = ones(n_controls, (horizon_disc - 1)); % size depends from the number of controls

    [~, usz] = size(u);
    u = [u(:, 2:end), zeros(n_controls, (horizon_disc - 1) - (usz - 1))];
    u=u*0;

    state_array = zeros(n_states, horizon_disc);
    state_array(:, 1) = istate;
    time_array = zeros(1, horizon_disc);
    time_array(1, 1) = it;
    t = it; %inital t

    for elem = 1:horizon_disc - 1 %compute the forward dynamics to define the error
        %dy = ForwardDynamics(state_array(:, elem), u(elem));
        %state_array(:, elem + 1) = euler_integration_fun(state_array(:, elem), dy, dt);
        state_array(:, elem + 1) = dynamics_rk4(state_array(:, elem), u(elem), dt);
        t = t + dt;
        time_array(1, elem + 1) = t;
    end
    

    tiledlayout(2, 1);
    nexttile
    plot(time_array, [state_array(1, :); state_array(3, :)])
    legend("x", "phi")
    nexttile
    plot(time_array, [state_array(2, :); state_array(4, :)])
    legend("dx", "dphi")
    pause(1e-4)

    %cost = fullcost(state_array, u, state_d, Q, R);

    while (max(abs(d(:, :))) > 1e-2)

        S(:, horizon_disc * 4 - 3:horizon_disc * 4) = Q * 1000;
        s(:, horizon_disc) = Q * 1000 * (state_array(:, horizon_disc) - state_d(:, horizon_disc));

        A_ = zeros(n_states, (horizon_disc - 1) * n_states);
        B_ = zeros(n_states, horizon_disc - 1);

        for n = (horizon_disc - 1):-1:1
            N = (4 * n - 3):(4 * n);

            %derivatives
            q = Q * (state_array(:, n) - state_d(:, n));
            r = R * u(n);
            [A_(:, N), B_(:, n)] = linearization_discretization(u(n), state_array(:, n),1);

            A = A_(:, N);
            B = B_(:, n);

            gx = q + A.' * s(:, n + 1);
            gu = r + B.' * s(:, n + 1);

            Gxx = Q + A.' * S(:, N + 4) * A;
            Guu = R + B.' * S(:, N + 4) * B; %H

            Gxu = A.' * S(:, N + 4) * B;
            Gux = B.' * S(:, N + 4) * A;

            d(:, n) = pinv(Guu) * gu; %d
            L(:, N) = pinv(Guu) * Gux; %K

            S(:, N) = Gxx + L(:, N).' * Guu * L(:, N) - Gxu * L(:, N) - L(:, N).' * Gux;
            s(:, n) = gx - L(:, N).' * gu  - Gxu * d(:, n);

            %s(:, n) = gx - L(:, N).' * gu - L(:, N).' * Guu * d(:, n) - Gxu * d(:, n);

            % dJ = dJ + gu' * d(:, n);

        end
        if any(isnan(L)) || any(any(isnan(d)))
            disp("L,l is nan")
            pause
        end

        new_state_array = state_array;
        new_u = u;

        %forward iteration
        Forward=false;
        while ~Forward
        for n = 1:horizon_disc - 1
            N = (4 * n - 3):(4 * n); %Prendere elemento da n*4-3 a elemento n*4

            new_u(n) = u(n) - alpha*d(:, n) - L(:, N) * (new_state_array(:, n) - state_array(:, n));

            %dy = ForwardDynamics(new_state_array(:, n), new_u(n));
            %new_state_array(:, n + 1) = euler_integration_fun(new_state_array(:, n), dy, dt);
            new_state_array(:, n + 1) = dynamics_rk4(new_state_array(:, n), new_u(n), dt);
        end
        if any(isnan(new_u)) || any(any(isnan(new_state_array))) || any(any(new_state_array>1e5))
            alpha=alpha/2;
        else
            Forward=true;
            alpha=0.5;
        end
        end

        %cost = fullcost(new_state_array, new_u, state_d, Q, R);

        state_array = new_state_array;
        u = new_u;

        tiledlayout(2, 1);
        nexttile
        plot(time_array, [state_array(1, :); state_array(3, :)])
        legend("x", "phi")
        nexttile
        plot(time_array, [state_array(2, :); state_array(4, :)])
        legend("dx", "dphi")
        pause(0.1)
        max(abs(d(:, :)))
       
    end

end

function c = stage_cost(x, u, x_d, Q, R)
    c = 0.5 * (x - x_d).' * Q * (x - x_d) + 0.5 * u.' * R * u;
end

function c = terminal_cost(x, x_d, Q)
    Qn = Q * 1000;
    c = 0.5 * (x - x_d).' * Qn * (x - x_d);
end

function J = fullcost(xtraj, utraj, x_d, Q, R)
    J = 0;
    [~, hor] = size(xtraj);

    for k = 1:hor - 1
        J = J + stage_cost(xtraj(:, k), utraj(k), x_d(:, k), Q, R);
    end

    k = hor;
    J = J + terminal_cost(xtraj(:, k), x_d(:, k), Q);
end

function state = dynamics_rk4(state, u, dt)
    f1 = ForwardDynamics(state, u);
    f2 = ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end
