function u = iLQR_function(istate, state_d, it, u)
    %global u %use global u to remember last optimized values
    dt = 0.01; %delta t for integration
    [n_states, sz] = size(state_d);
    n_controls = 1;


    %*define Q,R,Qn matrix
    Q = eye(n_states) * 0.01;
    Q(1, 1) = 5;
    Q(3, 3) = 1;
    R = 0.01;
    Qn = Q * 10000;

    iterations = 10000;
    bad_iterations=0;
    j_rm = 0.001;

    horizon = 3.00; %time S
    horizon_disc = floor(horizon / dt) + 1;
    defects_max = ones(n_states, 1) * 0.5; %difetti massimi per cui validare i controlli
    defects_max(1, 1) = 0.2;
    defects_max(3, 1) = 0.2;
    defects_max(2, 1) = 40;
    defects_max(4, 1) = 40;

    %* find real horizon
    if horizon_disc > (sz)
        horizon = (sz - 1) * dt;
        horizon_disc = (sz);
    end

    %*define desidered state in the horizon we have
    state_d = state_d(:, 1:horizon_disc);

    %*check and fix control lenght
    [~, usz] = size(u); 
    u = [u(:, 2:end), zeros(n_controls, (horizon_disc - 1) - (usz - 1))];
    
    state_array = zeros(n_states, horizon_disc);
    state_array(:, 1) = istate;
    new_u = u;
    new_state_array=state_array;

    defects = state_array(:, :) * 0;

    time_array = it:dt:(horizon + it);

    % !finished initialization
    
    for elem = 1:(horizon_disc - 1) %compute the forward dynamics to define the defects
        state_array(:, elem + 1)=dynamics_euler(state_array(:,elem),u(:,elem),dt);
    end

    J = cost(state_array, state_d, new_u, Q, R, Qn);
    new_J = J;
    disp("J")
    disp(J)

    % defects=state_array(:,:) - state_d(:,:);
    

    %start the optimizing iterations
    for iteration = 1:iterations - 1
        [L, l, A_, B_] = backward(n_states, horizon_disc, defects, state_array, state_d, u, Q, R, Qn);

        [new_state_array, new_u] = forward_linear_shoot(istate, horizon_disc, u, dt, L, l, state_array, defects,A_,B_,new_J);
        
        plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defects"],state_d,[1,3;2,4;5,NaN])
        pause

        % if any(isnan(new_u)) || any(any(new_state_array > [1e8; 1e8; 1e8; 1e8]))
        %     old_beta = beta;
        %     beta = beta / 2; %or beta=0
        %     alpha = alpha / 2;
        % else
        %     %alpha=old_alpha;
        %     solution_found = true;
        % end

        %[new_state_array,new_u] = forward_shoot(new_state_array(:, 1),horizon_disc,new_u,dt,L,l,state_array);

        [new_state_array, defects, new_u] = forward_multi_shoot(new_state_array(:, 1), horizon_disc, new_state_array, state_d, new_u, dt, L, l, state_array);

        %plot before exit to understand what is happening
        plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defects"],state_d,[1,3;2,4;5,NaN])

        %calculate new costs
        old_new_J = new_J;
        new_J = cost(new_state_array, state_d, new_u, Q, R, Qn);

        relative = abs(new_J - J) / new_J;
        if ~any(isnan(new_J)) && J > new_J
            bad_iterations=0;
        else
            bad_iterations=bad_iterations+1;
        end
        T = table(new_J, old_new_J, J, relative,bad_iterations, 'VariableNames', {'new_cost', 'prev_new_cost', 'prev_min_cost', 'relative','last bad'});
        disp(T)

        if all(~isnan(new_J)) && (J > new_J || relative < j_rm) %&& J ~= new_J
            state_array = new_state_array;
            J = new_J;
            u = new_u;
       
        end
        if (((relative < j_rm)) && all((sum(defects)) < defects_max))

            plot_xu(state_array, u, time_array, ["x", "dx", "phi", "dphi"],state_d,[1,3;2,4])


            return
        end
        pause

    end

end %function

%% Dynamics

function state = dynamics_rk4(state, u, dt)
    f1 = ForwardDynamics(state, u);
    f2 = ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end %function

function state = dynamics_euler(state, u, dt)
    state = euler_integration_fun(state, ForwardDynamics(state, u), dt);
end %function

%% Cost Function
function J = cost(state_array, state_d, u, Q, R, Qn)
    [~, horizon_disc] = size(state_array);
    J = 0;

    for i = 1:horizon_disc - 1
        J = J + (state_array(:, i) - state_d(:, i))' * Q * (state_array(:, i) - state_d(:, i)) + u(i)' * R * u(i);
    end

    i = horizon_disc;
    J = J + (state_array(:, i) - state_d(:, i))' * Qn * (state_array(:, i) - state_d(:, i));
end %function

%% Backward
function [L, l, A_, B_] = backward(n_states, horizon_disc, defects, x, xd, u, Q, R, Qn)
    S = zeros(n_states, horizon_disc * n_states);
    s = zeros(n_states, horizon_disc); %deep horizon+1 and hight is n_states
    L = zeros(1, horizon_disc - 1 * n_states); %size depends both from the number of controls and states
    l = zeros(1, horizon_disc - 1); % size depends from the number of controls
    %Fill the S and s matrix
    s(:, horizon_disc) = Qn * (x(:, horizon_disc) - xd(:, horizon_disc));
    S(:, horizon_disc * 4 - 3:horizon_disc * 4) = Qn;
    A_ = zeros(n_states, (horizon_disc - 1) * n_states);
    B_ = zeros(n_states, horizon_disc - 1);

    %back
    for n = (horizon_disc - 1):-1:1
        N = (4 * n - 3):(4 * n);
        [A_(:, N), B_(:, n)] = linearization_discretization(u(:, n), x(:, n), 1); %maybe discrete
        A = A_(:, N);
        B = B_(:, n);

        %compute r,h,G,H to simplify S and s computations
        P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
        r = R * u(:, n); % the derivative of uRu
        q = Q * (x(:, n) - xd(:, n));
        h = r + B.' * (s(:, n + 1) + S(:, N + 4) * (defects(:, n))); %gu + S(:, N + 4) * defects(:, n+1)
        G = P + B.' * S(:, N + 4) * A; %Gux
        H = R + B.' * S(:, N + 4) * B; %Guu
        L(:, N) = -pinv(H) * G; %K
        l(:, n) = -pinv(H) * h; %d

        Gxx = Q + A.' * S(:, N + 4) * A;
        gx = q + A' * (s(:, n + 1) + S(:, N + 4) * (defects(:, n)));

        % Gxu = G';
        % Gux = G;

        S(:, N) = Gxx - L(:, N).' * H * L(:, N);
        s(:, n) = gx + G.' * l(:, n) + L(:, N).' * (h + H * l(:, n));

    end

end %function

%% Forward Shoot
function [x,u] = forward_shoot(ix, horizon_disc, u, dt, L, l, x_old)

    x = repmat(ix, 1, horizon_disc);
    for n = 1:horizon_disc - 1
        N = (4 * n - 3):(4 * n);
        u(:, n) = u(:, n) + l(:, n) + L(:, N) * (x(:, n) - x_old(:, n));

        x(:, n + 1) = dynamics_rk4(x(:, n), u(n), dt);
    end
end %function

%% Forward Linear Shoot
function [x, u] = forward_linear_shoot(ix, horizon_disc, u, dt, L, l, x_old, defects,A_,B_,J)
    x = x_old * 0;
    x(:, 1) = ix; %initial state

    for n = 1:horizon_disc - 1
        N = (4 * n - 3):(4 * n);
        u(:, n) = u(:, n) + l(:, n) + L(:, N) * (x(:, n) - x_old(:, n));
        
        if J < 100000
            A = A_(:, N);
            B = B_(:, n);
            x(:, n + 1) = x_old(:, n + 1) ... % take the previous value
                + (A + B * L(:, N)) * (x(:, n) - x_old(:, n)) ... %like add the control
                +B * l(:, n) + (defects(:, n));
        else
            x(:, n + 1) = dynamics_rk4(x(:, n), u(:, n), dt);
        end
    end

end %function

%% Forward Multi Shoot
function [x, defects, u] = forward_multi_shoot(ix, horizon_disc, x_approx, state_d, u, dt, L, l, state_array)
    ix = ix(:);

    if horizon_disc < 40
        pieces = 1;
    else
        pieces = 1;
    end

    [~, len] = size(x_approx);
    len = floor(len / pieces);

    nn = zeros(4, pieces);
    ns = nn;

    for i = 1:pieces
        nn(:, i) = x_approx(:, len * (i - 1) + 1);
    end

    statess = zeros(length(ix), horizon_disc, pieces);
    uss = zeros(length(u(:, 1)), horizon_disc - 1, pieces);

    for i = 1:pieces

        stato = zeros(length(ix), horizon_disc);
        us = zeros(length(u(:, 1)), horizon_disc - 1);

        t = len * (i - 1);
        t = t + (t == 0); %skip if t==0
        inizio = t;

        stato(:, inizio) = nn(:, i);

        if i == pieces %if is last piece do till the end
            fine = horizon_disc - 1;
        else
            fine = len * (i) - 1;
        end

        for elem = inizio:fine
            n = t;
            N = (4 * n - 3):(4 * n);
            us(:, n) = u(:, n) + l(:, n) + L(:, N) * (stato(:, t) - state_array(:, n));
            stato(:, t + 1) = dynamics_rk4(stato(:, t), u(elem), dt);
            t = t + 1;
        end

        ns(:, i) = stato(:, t);

        if fine ~= horizon_disc - 1
            stato(:, t) = stato(:, t) * 0;
        end

        uss(:, :, i) = us(:, :);
        statess(:, :, i) = stato(:, :);
    end

    x = zeros(length(ix), horizon_disc);
    u = zeros(length(u(:, 1)), horizon_disc - 1);

    for i = 1:pieces
        x(:, :) = x(:, :) + statess(:, :, i);
        u(:, :) = u(:, :) + uss(:, :, i);
    end

    defects = zeros(length(ix), horizon_disc);

    for i = 1:1:pieces - 1
        defects(:, len * i) = ns(:, i) - nn(:, i + 1);
    end

    %defects(:, end) = ns(:, end) - state_d(:, end);
    defects = circshift(defects, -1, 2);
end %function

%% Plot

%plot_xu(x, u, time_array, names,xd,order)
function plot_xu(x, u, time_array, names,xd,order)
    if nargin < 3
        error("missing entries")
    end
     [h,~]=size(x);
    if nargin < 4
        names = repmat("", 1, h);
    end
    desired=true;
    if nargin < 5
        desired=false;
    end
    if nargin < 6

        order=(1:h)';
    end
    [h,l]=size(order);

    tiledlayout(h + 1, 1);

    for i = 1:h
        nexttile
        j=1;
        while j<=l && ~isnan(order(i,j)) 
            idx=order(i,j);
            plot(time_array, x(idx, :))
            legend(names(idx))
            if desired && length(xd(:,1))>=idx
                hold on
                plot(time_array, xd(idx, :),"--")
                legend(names(idx)+"_d")
            end 
            j=j+1;
        end
    end

    nexttile
    plot(time_array(1:end - 1), u)
    legend("controls")
    pause(1e-4)
end %function
