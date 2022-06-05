function next_single_control = iLQR_function(istate, state_d, it)
    global u %use global u to remember last optimized values
    next_single_control = 0; %default next control if we have a nan value
    dt = 0.01; %delta t for integration
    [n_states, sz] = size(state_d);
    n_controls = 1;
    alpha=0.5;

    Q = eye(n_states) * 0.001;
    Q(1, 1) = 1;
    Q(3, 3) = 1;
    R = 0.0001;
    Qn = Q * 1000;

    iterations = 10000;
    j_rm = 0.1;

    horizon = 5.99; %time S
    horizon_disc = floor(horizon / dt) + 1;
    defects_max = ones(n_states, 1) * 0.5; %difetti massimi per cui validare i controlli
    defects_max(1, 1) = 0.2;
    defects_max(3, 1) = 0.2;
    defects_max(2, 1) = 40;
    defects_max(4, 1) = 40;

    if horizon_disc > (sz)
        horizon = (sz - 1) * dt;
        horizon_disc = (sz);
    end

    state_d = state_d(:, 1:horizon_disc);

    S = zeros(n_states, horizon_disc * n_states);
    s = zeros(n_states, horizon_disc); %deep horizon+1 and hight is n_states
    L = zeros(1, horizon_disc - 1 * n_states); %size depends both from the number of controls and states
    l = zeros(1, horizon_disc - 1); % size depends from the number of controls

    if exist("u", "var")
        [~, usz] = size(u);
        u = [u(:, 2:end), zeros(n_controls, horizon_disc - usz)];
    else
        u = ones(n_controls, horizon_disc - 1) * (0);
    end

    new_u = u;

    state_array = zeros(n_states, horizon_disc);
    state_array(:, 1) = istate;
    time_array = it:dt:(horizon + it);

    for elem = 1:(horizon_disc - 1) %compute the forward dynamics to define the defects
        dy = ForwardDynamics(state_array(:, elem), u(elem));
        state_array(:, elem + 1) = euler_integration_fun(state_array(:, elem), dy, dt);
    end

    J = cost(state_array, state_d, new_u, Q, R, Qn);
    disp("J")
    disp(J)

    %cut defects in multiple pieces
    [~, len] = size(state_array);
    len = floor(len / 2);
    % defects=state_array(:,:) - state_d(:,:);
    defects = state_array(:, :) * 0;
    defects(:, end) = state_array(:, end) - state_d(:, end);
    defects = circshift(defects, -1, 2);
    
    %start the optimizing iterations
    for iteration = 1:iterations - 1

        %Fill the S and s matrix
        s(:, horizon_disc) = Qn * (defects(:, horizon_disc-1));
        S(:, horizon_disc * 4 - 3:horizon_disc * 4) = Qn;
        A_ = zeros(n_states, (horizon_disc - 1) * n_states);
        B_ = zeros(n_states, horizon_disc - 1);

        %back
        for n = (horizon_disc - 1):-1:1
            N = (4 * n - 3):(4 * n);
            [A_(:, N), B_(:, n)] = linearization_discretization(u(:, n), state_array(:, n), 1); %maybe discrete
            A = A_(:, N);
            B = B_(:, n);

            %compute r,h,G,H to simplify S and s computations
            P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
            r = R * u(:, n); % the derivative of uRu
            q = Q * (state_array(:, n) - state_d(:, n));
            h = r + B.' * (s(:, n + 1) + S(:, N + 4) * (defects(:, n))); %gu + S(:, N + 4) * defects(:, n+1)
            G = P + B.' * S(:, N + 4) * A; %Gux
            H = R + B.' * S(:, N + 4) * B; %Guu
            L(:, N) = -pinv(H) * G; %K
            l(:, n) = -pinv(H) * h; %d

            Gxx = Q + A.' * S(:, N + 4) * A;
            gx = q + A' * (s(:, n + 1) + S(:, N + 4) * (defects(:, n)));

            Gxu = G';
            Gux = G;
            S(:, N) = Gxx - L(:, N).' * H * L(:, N); %- Gxu * L(:, N) - L(:, N).' * Gux; %forse se so sbagliati con il meno
            s(:, n) = gx + G.' * l(:, n) + L(:, N).' * (h + H * l(:, n));

        end
        [L,l,A_,B_]=backward(horizon_disc,Q,R,Qn)

        new_state_array = zeros(n_states, horizon_disc);
        new_state_array(:, 1) = istate; %forward iteration
        new_u = zeros(n_controls, horizon_disc - 1);
        solution_found = false;
        alpha=1;
        while ~solution_found
            % disp("enter")
            for n = 1:horizon_disc - 1
                N = (4 * n - 3):(4 * n);
                new_u(:, n) = u(:, n) + alpha * l(:, n) +  L(:, N) * (new_state_array(:, n) - state_array(:, n));
                % A = A_(:, N);
                % B = B_(:, n);
                % new_state_array(:,n+1) =  state_array(:, n + 1) ... % take the previous value
                %     + (A + B * L(:, N)) * (new_state_array(:, n) - state_array(:, n)) ... %like add the control
                %     +B*alpha* l(:, n) + (defects(:,n));

                dy = ForwardDynamics(new_state_array(:, n), new_u(:, n));
                new_state_array(:, n + 1) = euler_integration_fun(new_state_array(:, n), dy, dt);

                %  state_array(:,elem+1)=dynamics_rk4(state_array(:, elem),new_u(elem),dt);

                % if isnan(new_state_array(:, n + 1))
                %     new_u(:, n) = NaN;
                %     disp("break")
                %     break
                % end
            end
            new_J=cost(new_state_array, state_d, new_u, Q, R, Qn);
            if new_J>J || isnan(new_J)
                alpha = alpha / 2;
                % disp(alpha)
            if alpha==0
                alpha=1;
            end
            else
                solution_found = true;
            end

        end

        state_array = zeros(n_states, horizon_disc);

        state_array(:, 1)=istate;
        for elem = 1:horizon_disc - 1
            state_array(:,elem+1)=dynamics_rk4(state_array(:, elem),new_u(elem),dt);
        %     %compute the forward dynamics to define the defects
            % dy = ForwardDynamics(state_array(:, elem), new_u(elem));
            % state_array(:, elem + 1) = euler_integration_fun(state_array(:, elem), dy, dt);
        end
        
        if horizon_disc<8
            pieces=1;
        else
            pieces=4;
        end

        [~, len] = size(state_array);
        len = floor(len / pieces);

        nn=zeros(4,pieces);
        ns=nn;
        for i=1:pieces
            nn(:,i)=state_array(:,len*(i-1)+1);
        end
        
        parfor i=1:pieces
            stato=nn(:,i);
            if i==pieces
                fine=horizon_disc-1;
            else
                fine=len*(i)-1;
            end
            for elem = len*(i-1):fine
                if elem==0 
                    continue 
                end
            %   stato=dynamics_rk4(stato,new_u(elem),dt);
            %   compute the forward dynamics to define the defects
                dy = ForwardDynamics(stato, new_u(elem));
                stato = euler_integration_fun(stato, dy, dt);
            end
            ns(:,i)=stato;
        end
        

        for i=1:pieces-1
            defects(:, len*i) = ns(:,i)- nn(:, i+1) ;
        end
        defects(:, end) = ns(:,end) - state_d(:, end);
        
%         defects(:, len) = state_array(:, len) - state_d(:, len);
%         defects(:, end) = state_array(:, end) - state_d(:, end);

        defects = circshift(defects, -1, 2);

        new_J = cost(state_array, state_d, new_u, Q, R, Qn);

        relative = abs(new_J - J) / new_J;
        disp('new_cost')
        disp(new_J)
        if J>new_J
            J=new_J;
            u = new_u;
            alpha=alpha*2;
        end
        next_single_control = u(1);
        if (((relative < j_rm)) && all((abs(state_array(:, horizon_disc) - state_d(:, horizon_disc))) < defects_max))

            %plot before exit to understand what is happening
            tiledlayout(2, 1);
            nexttile
            plot(time_array, [state_array(1, :); state_array(3, :)])
            legend("x", "phi")
            nexttile
            plot(time_array, [state_array(2, :); state_array(4, :)])
            legend("dx", "dphi")
             pause(0.01)
           
            return
        end
        
    end
end

function state = dynamics_rk4(state, u, dt)
    f1 = ForwardDynamics(state, u);
    f2 = ForwardDynamics(state + 0.5 * dt * f1, u);
    f3 = ForwardDynamics(state + 0.5 * dt * f2, u);
    f4 = ForwardDynamics(state + dt * f3, u);
    state = state + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end

function J = cost(state_array, state_d, u, Q, R, Qn)
    [~, horizon_disc] = size(state_array);
    J = 0;

    for i = 1:horizon_disc - 1
        J = J + (state_array(:, i) - state_d(:, i))' * Q * (state_array(:, i) - state_d(:, i)) + u(i)' * R * u(i);
    end

    i = horizon_disc;
    J = J + (state_array(:, i) - state_d(:, i))' * Qn * (state_array(:, i) - state_d(:, i));
end
