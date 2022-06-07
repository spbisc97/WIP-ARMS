function u = iLQR_function(istate, state_d, it,u)
    %global u %use global u to remember last optimized values
    next_single_control = 0; %default next control if we have a nan value
    dt = 0.01; %delta t for integration
    [n_states, sz] = size(state_d);
    n_controls = 1;
    alpha=0.5;
    beta=1;


    Q = eye(n_states) * 0.1;
    Q(1, 1) = 5;
    Q(3, 3) = 1;
    R = 0.0001;
    Qn = Q * 100;

    iterations = 10000;
    j_rm = 0.001;

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
    new_J=J;
    disp("J")
    disp(J)

    
    % defects=state_array(:,:) - state_d(:,:);
    defects = state_array(:, :) * 0;
    defects(:, end) = state_array(:, end) - state_d(:, end);
    defects = circshift(defects, -1, 2);
    
    %start the optimizing iterations
    for iteration = 1:iterations - 1

        [L,l,A_,B_]=backward(n_states,horizon_disc,defects,state_array,state_d,u,Q,R,Qn);


        new_state_array = zeros(n_states, horizon_disc);
        new_state_array(:, 1) = istate; %forward iteration
        new_u = zeros(n_controls, horizon_disc - 1);
        solution_found = false;
        old_beta=beta;
        while ~solution_found
            % disp("enter")
            for n = 1:horizon_disc - 1
                N = (4 * n - 3):(4 * n);
                new_u(:, n) = u(:, n) + alpha * l(:, n) +  beta *L(:, N) * (new_state_array(:, n) - state_array(:, n));
                if J<0
                A = A_(:, N);
                B = B_(:, n);
                new_state_array(:,n+1) = state_array(:, n + 1) ... % take the previous value
                    + (A + B * beta * L(:, N)) * (new_state_array(:, n)-state_array(:, n)) ... %like add the control
                    +B* alpha* l(:, n) + (defects(:,n));
                else
                dy = ForwardDynamics(new_state_array(:, n), new_u(:, n));
                new_state_array(:, n + 1) = euler_integration_fun(new_state_array(:, n), dy, dt);
                end
                %  state_array(:,elem+1)=dynamics_rk4(state_array(:, elem),new_u(elem),dt);

                % if isnan(new_state_array(:, n + 1))
                %     new_u(:, n) = NaN;
                %     disp("break")
                %     break
                % end
            end
            if  any(isnan(new_u)) || any(any(new_state_array>[1e8;1e8;1e8;1e8]))
                old_beta=beta;
                beta=beta/2; %or beta=0
                alpha=alpha/2;
            else
                beta=old_beta;
                %alpha=old_alpha;
                solution_found = true;
            end

        end
        %[new_state_array,defects,new_u] = forward_shoot(new_state_array(:, 1),horizon_disc,state_d,new_u,dt,L,l,state_array,alpha);

        [new_state_array,defects,new_u]=forward_multi_shoot(new_state_array(:, 1),horizon_disc,new_state_array,state_d,new_u,dt,L,l,state_array,alpha);

        %plot before exit to understand what is happening
            tiledlayout(3, 1);
            nexttile
            plot(time_array,[state_d(1,:);state_d(3,:);new_state_array(1, :); new_state_array(3, :)])
            legend("xd", "phid","x", "phi")
            nexttile
            plot(time_array, [new_state_array(2, :); new_state_array(4, :)])
            legend("dx", "dphi")
            nexttile
            plot(time_array(1:end-1),new_u)
            pause(0.01)

            old_new_J=new_J;
        new_J = cost(new_state_array, state_d, new_u, Q, R, Qn);
        if floor(old_new_J*1e1)==floor(new_J*1e1)
            beta=beta/3;
            disp("beta split")
            disp(beta)
        end

        relative = abs(new_J - J) / new_J;
        T=table(new_J,old_new_J,J,relative,'VariableNames',{'new_cost','prev_new_cost','prev_min_cost','relative'});
        % disp(['new_cost ',' min_cost'])
        % disp([new_J,J])
        % disp("relative")
        % disp(relative)
        disp(T)
        if all(~isnan(new_J))&&(J>new_J || relative < j_rm )&& J~=new_J
            state_array=new_state_array;
            J=new_J;
            u = new_u;
            alpha=alpha*2;
        else
            alpha=alpha/2;
        end
        next_single_control = u(1);
        if (((relative < j_rm)) && J<1e3 && all((abs(state_array(:, horizon_disc) - state_d(:, horizon_disc))) < defects_max))

            %plot before exit to understand what is happening
            tiledlayout(3, 1);
            nexttile
            plot(time_array, [ state_d(1,:);state_d(3,:);state_array(1, :); state_array(3, :)])
            legend("xd", "phid","x", "phi")
            nexttile
            plot(time_array, [state_array(2, :); state_array(4, :)])
            legend("dx", "dphi")
            nexttile
            plot(time_array(1:end-1),u)
            pause(0.01)
            return
        end
        
    end
end




%
%
%Auxiliary Functions

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




function [L,l,A_,B_] =backward(n_states,horizon_disc,defects,state_array,state_d,u,Q,R,Qn)
    S = zeros(n_states, horizon_disc * n_states);
    s = zeros(n_states, horizon_disc); %deep horizon+1 and hight is n_states
    L = zeros(1, horizon_disc - 1 * n_states); %size depends both from the number of controls and states
    l = zeros(1, horizon_disc - 1); % size depends from the number of controls
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
        q = Q * (state_array(:,n)-state_d(:,n));
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
end

function [x,defects,u] = forward_shoot(ix,horizon_disc,state_d,u,dt,L,l,state_array,alpha)
   
    x=repmat(ix,1,horizon_disc);
    for elem = 1:horizon_disc - 1
        n=elem;
        N = (4 * n - 3):(4 * n);
        u(:, n) = u(:, n) + alpha * l(:, n) +  L(:, N) * (x(:, n) - state_array(:, n));
        x(:,elem+1)=dynamics_rk4(x(:, elem),u(elem),dt);
    %     %compute the forward dynamics to define the defects
        % dy = ForwardDynamics(state_array(:, elem), new_u(elem));
        % state_array(:, elem + 1) = euler_integration_fun(state_array(:, elem), dy, dt);
    end
    defects=zeros(length(ix),horizon_disc);

    %defects(:, end) = x(:,end) - state_d(:, end);
    defects = circshift(defects, -1, 2);

end


function [x,defects,u]= forward_multi_shoot(ix,horizon_disc,x_approx,state_d,u,dt,L,l,state_array,alpha)  
    ix=ix(:);
    if horizon_disc<40
        pieces=1;
    else
        pieces=8;
    end

    [~, len] = size(x_approx);
    len = floor(len / pieces);

    nn=zeros(4,pieces);
    ns=nn;
    for i=1:pieces
        nn(:,i)=x_approx(:,len*(i-1)+1);
    end
    statess=zeros(length(ix),horizon_disc,pieces);
    uss=zeros(length(u(:,1)),horizon_disc-1,pieces);

    parfor i=1:pieces

        stato=zeros(length(ix),horizon_disc);
        us=zeros(length(u(:,1)),horizon_disc-1)
        t=len*(i-1);
        t=t+(t==0);
        stato(:,t)=nn(:,i);
        % statess(:,i,t)=stato;
        if i==pieces %if is last piece do till the end
            fine=horizon_disc-1;
        else
            fine=len*(i)-1;
        end
        for elem = len*(i-1):fine
            if elem==0 
                continue 
            end

        n=t;
        N = (4 * n - 3):(4 * n);
        us(:, n) = u(:, n) + alpha * l(:, n) +  L(:, N) * (stato(:,t) - state_array(:, n));

            stato(:,t+1)=dynamics_rk4(stato(:,t),u(elem),dt);
        %   compute the forward dynamics to define the defects
%             dy = ForwardDynamics(stato(:,t), u(elem));
%             stato(:,t+1) = euler_integration_fun(stato(:,t), dy, dt);
            t=t+1;
        end
        ns(:,i)=stato(:,t);
        if fine~=horizon_disc-1
        stato(:,t)=stato(:,t)*0;
        end
        uss(:,:,i)=us(:,:);
        statess(:,:,i)=stato(:,:);
    end

    x=zeros(length(ix),horizon_disc);
    u=zeros(length(u(:,1)),horizon_disc-1);

    for i=1:pieces
        x(:,:)=x(:,:)+statess(:,:,i);
        u(:,:)=u(:,:)+uss(:,:,i);
    end
    defects=zeros(length(ix),horizon_disc);
    for i=1:pieces-1
        defects(:, len*i) = ns(:,i)-nn(:, i+1) ;
    end
    defects(:, end) = ns(:,end) - state_d(:, end);
%     tiledlayout(1,1)
%     nexttile
%     plot(1:horizon_disc,x(:,:))
% 
%     pause
    defects = circshift(defects, -1, 2);

end