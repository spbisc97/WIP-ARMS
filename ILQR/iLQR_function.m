function next_single_control = iLQR_function(istate, state_d, it)

    % [~, sz] = size(state_d) ;

    % state_d = reshape(state_d, 4, sz/4);

    % disp(state_d)
    global u
    next_single_control = 0;

    dt = 0.01;
    t = it;
    state = istate;
    %disp("ilqr")

    %pause

    [n_states, sz] = size(state_d);

    Q = eye(n_states) * 0.1;
    Q(1, 1) = 30;
    Q(3, 3) = 10;
    R = 0.0001;

    iterations = 10000;
    cost = 0;
    new_cost = 0;
    j_rm = 0.1;

    horizon = 0.3 ; %time S
    horizon_disc = floor(horizon / dt);
    defects_max = ones(n_states, 1) * 0.3;
    defects_max(1, 1) = 0.2;
    defects_max(2, 1) = 30;
    defects_max(4, 1) = 30;

    if horizon_disc > sz
        horizon = sz * dt;
        horizon_disc = sz;
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

    % myvals=whos;
    % for n = 1:length(myvals)
    %     if isnan(myvals(n).name)
    %       eval(myvals(n).name)
    %     end
    % end

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
    %defects(:,1:end-1)=0;
    %disp("defects")
    %disp(defects(:)')
    s = zeros(n_states, horizon_disc); %deep horizon+1 and hight is n_states
    s(:, horizon_disc) = Q * defects(:, horizon_disc); %TBC % check if it is Q*defect(horizon disc)

    %start the optimizing iterations
    for iteration = 1:iterations - 1
        %fill the s matrix , expecially the last element
        %last element = horizon_disc
        A_ = [];
        B_ = [];
        %compute the linearized dynamics
        for step = 1:horizon_disc-1
            Step = (4 * step - 3):(4 * step);
            [A_(:, Step), B_(:, step)] = linearization_discretization(u(step), state_array(:, step));
            %A_=[A_,A]
            %B_=[B_,B]
        end

        %backword iteration
        for step = 1:(horizon_disc-1) %horizon_disc-1 times
            n = horizon_disc - step;
            N = (4 * n - 3):(4 * n); %Prendere elemento da n*4-3 a elemento n*4
            A = A_(:, N);
            B = B_(:, n);

            %compute r,h,G,H to simplify S and s computations
            P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
            r = R * u(n); %should be the derivative of uRu
            h = r + B.' * (s(:, n + 1) + S(:, N + 4) * defects(:, n));
            G = P + B.' * S(:, N + 4) * A;
            H = R + B.' * S(:, N + 4) * B;
            % disp(["H"])
            % disp([H])
            % disp("G")
            % disp([G])
            % disp("h")
            % disp([h])
            %compute Values to use in the forward iterations
            L(:, N) = -pinv(H) * G;
            l(:, n) = -pinv(H) * h;
            %disp(["L","l","n"])
            %disp([L(n),l(n),n])

            %compute next S,s Value
            S(:, N) = Q + A' * S(:, N + 4) * A - L(:, N)' * H * L(:, N);
            q = Q * defects(:, n);
            s(:, n) = q + A' * (s(:, n + 1) + S(:, N + 4) * defects(:, n)) + G' * l(:, n) + L(:, N)' * (h + H * l(:, n));
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
        new_state_array=[istate];
        %forward iteration
        for n = 1:horizon_disc-1
            %compute delta_u and update the value
            N = (4 * n - 3):(4 * n); %Prendere elemento da n*4-3 a elemento n*4
            A = A_(:, N);
            B = B_(:, n);

            delta_u = l(:, n) + L(:, N) *(new_state_array(:,n)-state_array(:,n)); %defects(:, n);
            new_u(n) = u(n) + delta_u;
            new_state=state_array(:,n+1)+(A+B*L(:,N))*(new_state_array(:,n)-state_array(:,n))...
                                    +B*l(:,n)+defects(:,n+1);
            new_state_array=[new_state_array,new_state];
            if (isnan(new_u(n)))
                disp(iteration)
                disp("isnan")
                pause
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
        %defects(:,1:end-1)=0;
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
        tiledlayout(2,1);
        nexttile
        plot(time_array,[state_array(1,:);state_array(3,:)])
        legend("x", "phi")
        nexttile
        plot(time_array, [state_array(2,:); state_array(4,:)])
        legend("dx", "dphi")

        new_cost = 0;
        
        for i = 1:horizon_disc-1
            new_cost = new_cost + defects(:, i+1)' * Q * defects(:, i+1) + control_array(:, i)' * R * control_array(:, i);
        end
        relative=abs(new_cost - cost) / new_cost;
        disp('new_cost')
        disp(new_cost)
        disp('relative')
        disp(relative)

        pause(0.5)
        u = new_u;
        next_single_control = u(1);

        %check cost increments and return if solved

        if (((relative < j_rm)) && all((abs(defects(:,horizon_disc)'))'<defects_max))
            
            return
        end

        cost = new_cost;

        

    end

    % disp(u(1))
    
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
