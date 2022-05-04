function u_next = iLQR_function(istate, state_d, it)
    
    % [~, sz] = size(state_d) ;
    
    % state_d = reshape(state_d, 4, sz/4);
    
    % disp(state_d)
    global u
    u_next = 0;

    dt = 0.01;
    t = it;
    state = istate;
    %disp("ilqr")

    %pause

    [n_states, sz] = size(state_d);

    Q = eye(n_states) * 20;
    Q(1,1)=70;
    Q(3,3)=200;
    R = 0.01;

    iterations = 100;
    cost=0;
    new_cost=0;
    j_rm=0.002;
    defects_max=ones(n_states,1);
    defects_max(2,1)=20;
    defects_max(4,1)=20;
    horizon = 3; %time S
    horizon_disc = floor(horizon / dt);

    if horizon_disc > sz
        horizon = sz * dt;
        horizon_disc = sz;
    end

    S = repmat(Q, 1, horizon_disc + 1);
    L = zeros(1, horizon_disc * n_states); %size depends both from the number of controls and states
    l = zeros(1, horizon_disc); % size depends from the number of controls

    if exist("u", "var")
        [usz, ~] = size(u);
        u = [0; u];
        u = [u(3:end); zeros(horizon_disc - usz + 1, 1)];
        next_u=u;
    else
        u = ones(horizon_disc, 1) * (0);
    end

    state_array = [];
    control_array = [];
    time_array = [];

    % myvals=whos;
    % for n = 1:length(myvals)
    %     if isnan(myvals(n).name)
    %       eval(myvals(n).name)
    %     end
    % end

    while t < it + horizon
        %compute the forward dynamics to define the defects
        elem = floor((t - it + dt) / dt);
        dy = ForwardDynamics(state, u(elem));
        state = euler_integration_fun(state, dy, dt);
        control_array = [control_array, u(elem)];
        time_array = [time_array, t];
        state_array = [state_array, state];
        t = t + dt;
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
    defects = state_array(:, 1:horizon_disc) - state_d(:, 1:horizon_disc);
    % defects(1:end-1)=0;
    %disp("defects")
    %disp(defects(:)')
    s = zeros(n_states, horizon_disc + 1); %deep horizon+1 and hight is n_states
    s(:, horizon_disc + 1) = Q * state_array(:, horizon_disc); %TBC

    %start the optimizing iterations
    for iteration = 1:iterations - 1
        %fill the s matrix , expecially the last element
        %last element = horizon_disc
        A_ = [];
        B_ = [];
        %compute the linearized dynamics
        for step = 1:horizon_disc
            Step=(4*(step-1)+1):(4*(step));
            [A_(:,Step), B_(:,step)] = linearization_discretization(u(step), state_array(:,step));
        end

        %backword iteration
        for step = 1:horizon_disc %horizon_disc-1 times
            n = horizon_disc - step + 1; 
            N = (4*n-3):(4*n); %Prendere elemento da n*4-3 a elemento n*4
            A = A_(:,N);
            B = B_(:,n);

            %compute r,h,G,H to simplify S and s computations
            P = 0;%repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
            r = R * u(n); %should be the derivative of uRu
            h = r + B.' * (s(:,n + 1) + S(:,N + 4) * defects(:, n));
            G = P + B.' *  S(:,N + 4) * A;
            H = R + B.' *  S(:,N + 4) * B;
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
            S(:,N) = Q + A' * S(:,N + 4) * A - L(:, N)' * H * L(:, N);
            q = Q * defects(:, n);
            s(:,n) = q + A' * (s( :,n + 1) + S(:,N + 4) * defects(:, n)) + G' * l(:, n) + L(:, N)' * (h + H * l(:, n));
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

        %forward iteration
        for n = 1:horizon_disc
            %compute delta_u and update the value
            delta_u = l(:, n) + L(:, N) * defects(:, n);
            next_u(n) = u(n) + delta_u;
            if (isnan(next_u(n)))
                exit()
            end
        end

        t = it;
        state = istate;
        state_array = [];
        control_array = [];
        time_array = [];

        while t < it + horizon
            %compute the forward dynamics to define the defects
            elem = floor((t - it + dt) / dt);
            dy = ForwardDynamics(state, next_u(elem));
            state = euler_integration_fun(state, dy, dt);
            control_array = [control_array, u(elem)];
            time_array = [time_array, t];
            state_array = [state_array, state];
            t = t + dt;
        end

        %defects=distance between state array and desired state
        %(traj and desired trajectory)
        defects = state_array(:, 1:horizon_disc) - state_d(:, 1:horizon_disc);
        %defects(1:end-1)=0;
        disp("istate")
        disp(istate)
        disp("state array")
        disp(state_array(:,:))
        disp("state_d traj")
        disp(state_d(:,1:horizon_disc))
        disp("defects")
        disp(defects(:,:))
        disp("control")
        disp(u')
        plot(time_array, state_array)
        legend("x","dx","phi","dphi")

        new_cost=0;
        for i=1:horizon_disc
            new_cost=new_cost+defects(:,i)'*Q*defects(:,i)+control_array(:,i)'*R*control_array(:,i);
        end
        disp("new_cost")
        disp(new_cost)
        pause(0.1)
            %check cost increments and return if solved
            
        if((abs(new_cost-cost)/new_cost < j_rm) && sum(sum(abs(defects'))'<defects_max))
            return
        end
        
        cost=new_cost;

        




        u=u_next;

    end

    % disp("next control")
    % disp(u(1))
    u_next = u(1);
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
