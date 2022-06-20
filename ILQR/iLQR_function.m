function u = iLQR_function(istate, state_d, it, u)
    %global u %use global u to remember last optimized values
    dt = 0.01; %delta t for integration
    [n_states, sz] = size(state_d);
    n_controls = 1;    


    %*define Q,R,Qn matrix
    Q = eye(n_states) * 0.01;
    Q(1, 1) = 1;
    Q(3, 3) = 1;
    R = 0.001;
    Qn = Q*1;

    iterations = 100000;
    plot_step=10;
    bad_iterations=0;
    j_rm = 0.01;

    horizon = 5.7; %time S
    horizon_disc = floor(horizon / dt) + 1;
    defects_max = ones(n_states, 1) ; %difetti massimi per cui validare i controlli
        defects_max(1, 1) = 0.01;
        defects_max(3, 1) = 0.01;
        defects_max(2, 1) = 0.02;
        defects_max(4, 1) = 0.02;
    %defects_max=defects_max*2 

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
    % u=u*0; reset suggested controls
    state_array = ones(1, horizon_disc).*istate;
    new_u = u;
    alpha = 1;
    lmb=1;
    time_array = it:dt:(horizon + it);

    % !finished initialization

    %     for elem = 1:(horizon_disc - 1) %compute the forward dynamics to define the defects
    %         state_array(:, elem + 1)=dynamics_euler(state_array(:,elem),u(:,elem),dt);
    %     end


    L = zeros(1, horizon_disc  * n_states); %size depends both from the number of controls and states
    l = zeros(1, horizon_disc ); % size depends from the number of controls
    [state_array, defects, u ]= forward_multi_shoot(horizon_disc,state_array, u, dt,state_array,L,l);

    %plot_xu([state_array;defects], u, time_array, ["x", "dx", "phi", "dphi","defx","defdx","defphi","defdphi"],state_d,[1,3,nan,nan;2,4,nan,nan;5,6,7,8],"multi start")

    J = cost(state_array, state_d, new_u, Q, R, Qn)*1e3;
    new_J = J;
    disp("J")
    disp(J)


    %start the optimizing iterations
    for iteration = 1:iterations - 1
        [L, l, A_, B_] = backward(n_states, horizon_disc, defects, state_array, state_d, u, Q, R, Qn);

        [new_state_array, new_u] = forward_linear_shoot(horizon_disc,state_array, u, dt, lmb*L, alpha*l, defects,A_,B_,new_J);

        %plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defects"],state_d,[1,3;2,4;5,NaN],"linear")
                if mod(iteration,plot_step)==0

        plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defx","defdx","defphi","defdphi"],state_d,[1,3,nan,nan;2,4,nan,nan;5,6,7,8],"linear")
                end

        %[new_state_array,new_u] = forward_shoot(new_state_array(:, 1),horizon_disc,new_u,dt,L,l,state_array);

        [new_state_array, new_defects, new_u] = forward_multi_shoot(horizon_disc, new_state_array, new_u, dt,state_array, lmb*L, alpha*l);

        %plot before exit to understand what is happening
        %plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defects"],state_d,[1,3;2,4;5,NaN],"multi")
        if mod(iteration,plot_step)==0
        plot_xu([new_state_array;new_defects], new_u, time_array, ["x", "dx", "phi", "dphi","defx","defdx","defphi","defdphi"],state_d,[1,3,nan,nan;2,4,nan,nan;5,6,7,8],"multi")
        end
        %pause
        %calculate new costs
        old_new_J = new_J;
        new_J = cost(new_state_array, state_d, new_u, Q, R, Qn);

        relative = abs(new_J - J) / new_J;
        if ~isnan(new_J) && ~isinf(new_J)   && (J > new_J)
            bad_iterations=0;
        else
            bad_iterations=bad_iterations+1;
        end
        T = table(new_J, old_new_J, J,max(sum(abs(new_defects),2)), relative,bad_iterations,alpha, 'VariableNames', {'new_cost', 'prev_new_cost', 'prev_min_cost', 'max(sum(def))','relative','last bad','alpha'});

        disp(T)



        if bad_iterations==0
            state_array = new_state_array;
            J = new_J;
            u = new_u;
            defects=new_defects;
            alpha=1;
            lmb=1;
        else
            if ~isnan(new_J) && ~isinf(new_J)
                alpha=alpha/2;
                lmb=1;
                state_array = new_state_array;
                J = new_J;
                u = new_u;
                defects=new_defects;
            else
            alpha=alpha/2;
            lmb=1;
                continue
            end
        end


        if bad_iterations > 20
            relative=0;

        else
        end
        if (((relative < j_rm)) && all((sum(abs(new_defects),2)) < defects_max))

            plot_xu([state_array;defects], u, time_array, ["x", "dx", "phi", "dphi","defx","defdx","defphi","defdphi"],state_d,[1,3,nan,nan;2,4,nan,nan;5,6,7,8],"final")
            return
        end


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
    S(:, horizon_disc * n_states - (n_states-1):horizon_disc * n_states) = Qn;
    A_ = zeros(n_states, (horizon_disc - 1) * n_states);
    B_ = zeros(n_states, horizon_disc - 1);

    %back
    for n = (horizon_disc - 1):-1:1
        N = (n_states * n - (n_states-1)):(n_states * n);
        [A_(:, N), B_(:, n)] = linearization_discretization(u(:, n), x(:, n), 1);
        A = A_(:, N);
        B = B_(:, n);

        %compute r,h,G,H to simplify S and s computations
        P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
        r = R * u(:, n); % the derivative of uRu
        q = Q * (x(:, n) - xd(:, n));
        h = r + B.' * (s(:, n + 1) + S(:, N + n_states) * (defects(:, n))); %gu + S(:, N + 4) * defects(:, n+1)
        G = P + B.' * S(:, N + n_states) * A; %Gux
        H = R + B.' * S(:, N + n_states) * B; %Guu
        L(:, N) = -pinv(H) * G; %K
        l(:, n) = -pinv(H) * h; %d

        Gxx = Q + A.' * S(:, N + n_states) * A;
        gx = q + A' * (s(:, n + 1) + S(:, N + n_states) * (defects(:, n)));

        % Gxu = G';
        % Gux = G;

        S(:, N) = Gxx - L(:, N).' * H * L(:, N);
        s(:, n) = gx + G.' * l(:, n) + L(:, N).' * (h + H * l(:, n));

    end

end %function

%% Forward Shoot
function [x,u] = forward_shoot(ix, horizon_disc, u, dt, L, l, x_old)
    n_states=lenght(ix);
    x = repmat(ix, 1, horizon_disc);
    for n = 1:horizon_disc - 1
        N = (n_states * n - (n_states-1)):(n_states* n);
        u(:, n) = u(:, n) + l(:, n) + L(:, N) * (x(:, n) - x_old(:, n));

        x(:, n + 1) = dynamics_rk4(x(:, n), u(n), dt);
    end
end %function

%% Forward Linear Shoot
function [x, u] = forward_linear_shoot(horizon_disc,x_old, u, dt, L, l,  defects,A_,B_,J)
    x = x_old;
    [n_states,~]=size(x);


    %x(:, 1) = ix; %initial state

    for n = 1:horizon_disc - 1
        N = (n_states * n - (n_states-1)):(n_states * n);
        c=1;%c=all(defects(:, n)==0);
        %u(:, n) = u(:, n) + c*(l(:, n) + L(:, N) * (x(:, n) - x_old(:, n)));

        if J < 1e8
            A = A_(:, N);
            B = B_(:, n);
            x(:, n + 1) = x_old(:, n + 1)+(defects(:, n))... % take the previous value
                +c*((A + B * L(:, N)) * (x(:, n) - x_old(:, n)) ... %like add the control
                +B * l(:, n));

        else
            x(:, n + 1) = dynamics_euler(x(:, n), u(:, n), dt);
        end
    end

end %function

%% Forward Multi Shoot

function [x, defects, u] = forward_multi_shoot(horizon_disc,x, u, dt,x_old, L, l)


    %closedloop=0;

    n_states=length(x(:,1));
    n_controls=length(u(:,1));
    if horizon_disc < 3
        pieces = 1;
        k=1;
    else
        pieces = 16;
        %k=2.^(0:-0.1:(-pieces+1)/10);
        k=ones(1,pieces)*1;
    end

    [~, len] = size(x);
    len = floor(len / pieces);

    x_start = zeros(n_states, pieces);
    x_arrive = x_start;

    for i = 1:pieces
        x_start(:, i) = x(:, len * (i - 1) + 1);
    end

    statess = zeros(n_states, horizon_disc, pieces);
    uss = zeros(n_controls, horizon_disc - 1, pieces);



    parfor i = 1:pieces

        stato = zeros(n_states, horizon_disc);
        us = zeros(n_controls, horizon_disc - 1);

        t = len * (i - 1)+1;
        t = t + (t == 0); %skip if t==0
        inizio = t;

        stato(:, inizio) = x_start(:, i);

        if i == pieces %if is last piece do till the end
            fine = horizon_disc - 1;
        else
            fine = len * (i);
        end
        %len(i-1)+1:len(i):
        for elem = inizio:fine
            n = t;
            N = (n_states * n - (n_states-1)):(n_states * n);
            us(:, n) = u(:, n) +(k(i)*l(:, n)+  L(:, N) * (stato(:, t) - x_old(:, n)));
            stato(:, t + 1) = dynamics_rk4(stato(:, t), us(elem), dt);
            t = t + 1;
        end

        x_arrive(:, i) = stato(:, t);

        if fine ~= (horizon_disc - 1)
            stato(:, t) = stato(:, t) * 0;
        end

        uss(:, :, i) = us(:, :);
        statess(:, :, i) = stato(:, :);
    end

    x = zeros(n_states, horizon_disc);
    u = zeros(n_controls, horizon_disc - 1);

    for i = 1:pieces
        x(:, :) = x(:, :) + statess(:, :, i);
        u(:, :) = u(:, :) + uss(:, :, i);
    end

    defects = zeros(n_states, horizon_disc);

    for i = 1:1:pieces - 1
        defects(:, len*i+1) = (x_arrive(:, i) - x_start(:, i + 1));
    end
    %defects(:, end) = x_arrive(:, end) - state_d(:, end);
    defects = circshift(defects, -1, 2);
end %function

%% Plot


function plot_xu(x, u, time_array, names,xd,order,type)
    %plot_xu(x, u, time_array, names,xd,order)
    if type==""
        return
    end
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
    if nargin < 7
        type='figure';
    end
    [h,l]=size(order);

    tiledlayout(h + 1, 1);

    for i = 1:h
        nexttile
        j=1;
        lgd=[];
        while j<=l && ~isnan(order(i,j))
            idx=order(i,j);
            plot(time_array, x(idx, :))
            hold on
            lgd=[lgd,names(idx)]; %#ok
            if desired && length(xd(:,1))>=idx
                plot(time_array, xd(idx, :),LineStyle="--")
                hold on
                lgd=[lgd,names(idx)+"_d"];%#ok
            end
            title(type)

            j=j+1;
        end
        legend(lgd);

    end

    nexttile
    plot(time_array(1:end - 1), u)
    legend("controls")
    title(type)

    pause(1e-4)
    set(gcf,'Name',type,'NumberTitle','off')
end %function
