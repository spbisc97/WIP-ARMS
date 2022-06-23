classdef iLQR_GNMS
    properties
        Q double = 1
        R double = 0.1
        Qn double =1
        pieces  =16
        defects_max double =0.001
        horizon double = 10;
        horizon_disc=1001
        order=[]
        names=[]
        model
        plot_steps  = 1;
        plot_start  =true ;
        plot_end = true ;
        dt=0.01
    end
    methods
        %% Dynamics

        function state = dynamics_rk4(obj,state, u)
            f1 = obj.model.ForwardDynamics(state, u);
            f2 = obj.model.ForwardDynamics(state + 0.5 * obj.dt * f1, u);
            f3 = obj.model.ForwardDynamics(state + 0.5 * obj.dt * f2, u);
            f4 = obj.model.ForwardDynamics(state + obj.dt * f3, u);
            state = state + (obj.dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
        end %function

        function state = dynamics_euler(obj,state, u)
            state = obj.model.euler_integration_fun(state, ForwardDynamics(state, u), obj.dt);
        end %function

        function obj=set.horizon(obj,horizon)
            if horizon<0 
                error("horizon has to be positive")
            else
                obj.horizon=horizon;
                obj.horizon_disc=floor(obj.horizon/obj.dt)+1 ;%ok
            end
        end


        function obj=iLQR_GNMS(model,Q ,R ,Qn,pieces,defects_max)

            if nargin < 1
                error("Model is needed")
            end
            obj.model=model;
            if nargin > 2
                obj.Q=Q;
                obj.R=R;
            end
            obj.Qn=Q;
            if nargin > 3
                obj.Qn=Qn;
            end
            if nargin > 4
                obj.pieces=pieces;
            end
            if nargin > 5
                obj.defects_max=defects_max;
            end
        end
        function u = iLQR_function(obj,istate, state_d, it, u)
            %global u %use global u to remember last optimized values
            [n_states, sz] = size(state_d);
            [n_controls, ~] = size(u);

            iterations = 100000;
            bad_iterations = 0;
            j_rm = 0.0001;

            %* find real horizon
            if obj.horizon_disc > (sz)
                obj.horizon = (sz - 1) * obj.dt;
                obj.horizon_disc = (sz);
            end

            %*define desidered state in the horizon we have
            state_d = state_d(:, 1:obj.horizon_disc);

            %*check and fix control lenght
            [~, usz] = size(u);
            u = [u(:, 2:end), zeros(n_controls, (obj.horizon_disc - 1) - (usz - 1))];
            % u=u*0; reset suggested controls
            state_array = ones(1, obj.horizon_disc) .* istate;
            new_u = u;
            alpha = 1;
            lmb = 1;
            time_array = it:obj.dt:(obj.horizon + it);

            % !finished initialization

            %     for elem = 1:(horizon_disc - 1) %compute the forward dynamics to define the defects
            %         state_array(:, elem + 1)=dynamics_euler(state_array(:,elem),u(:,elem),obj.dt);
            %     end

            L = zeros(n_controls,n_states, obj.horizon_disc ); %size depends both from the number of controls and states
            l = zeros(n_controls, obj.horizon_disc); % size depends from the number of controls
            [state_array, defects, u] = forward_multi_shoot(obj, state_array, u, state_array, L, l);

            if obj.plot_start
                obj.plot_xu([state_array;defects], u, time_array, obj.names,state_d,obj.order,"multi start")
            end
            J = obj.cost(state_array, state_d, new_u) * 1e3;
            new_J = J;
            disp("J")
            disp(J)

            %start the optimizing iterations
            for iteration = 1:iterations - 1
                [L, l, A_, B_] = backward(obj,n_states,  defects, state_array, state_d, u);
                
                [new_state_array, new_u] = forward_linear_shoot(obj, state_array, u,lmb * L, alpha * l, defects, A_, B_, new_J);

                %     plot(time_array,state_array)
                %     legend()
                % disp("new_linear")
                % pause
                %plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defects"],state_d,[1,3;2,4;5,NaN],"linear")
                if mod(iteration, obj.plot_steps) == 0

                    %plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defx","defdx","defphi","defdphi"],state_d,[1,3,nan,nan;2,4,nan,nan;5,6,7,8],"linear")
                    obj.plot_xu([new_state_array; defects], new_u, time_array, obj.names, state_d, obj.order, "linear")
                end

                %[new_state_array,new_u] = forward_shoot(new_state_array(:, 1),horizon_disc,new_u,dt,L,l,state_array);

                [new_state_array, new_defects, new_u] = obj.forward_multi_shoot(new_state_array, new_u, state_array, lmb*L, alpha*l);

                %plot before exit to understand what is happening
                %plot_xu([new_state_array;defects], new_u, time_array, ["x", "dx", "phi", "dphi","defects"],state_d,[1,3;2,4;5,NaN],"multi")
                if mod(iteration, obj.plot_steps) == 0
                    %plot_xu([new_state_array;new_defects], new_u, time_array, ["x", "dx", "phi", "dphi","defx","defdx","defphi","defdphi"],state_d,[1,3,nan,nan;2,4,nan,nan;5,6,7,8],"multi")
                    obj.plot_xu([new_state_array; new_defects], new_u, time_array,obj.names, state_d,obj.order, "multi")

                end

                %pause
                %calculate new costs
                new_J = obj.cost(new_state_array, state_d, new_u);

                relative = abs(new_J - J) / new_J;

                if ~isnan(new_J) && ~isinf(new_J) && (J > new_J)
                    bad_iterations = 0;
                else
                    bad_iterations = bad_iterations + 1;
                end

                T = table(new_J, J, max(sum(abs(new_defects), 2)), relative, bad_iterations, alpha, 'VariableNames', {'new_cost', 'prev_cost', 'max(sum(def))', 'relative', 'last bad', 'alpha'});

                disp(T)

                if bad_iterations == 0
                    state_array = new_state_array;
                    J = new_J;
                    u = new_u;
                    defects = new_defects;
                    alpha = 1;
                    lmb = 1;
                else

                    if ~isnan(new_J) && ~isinf(new_J)
                        alpha = alpha / 2;
                        lmb = 1;
                        state_array = new_state_array;
                        J = new_J;
                        u = new_u;
                        defects = new_defects;
                    else
                        alpha = alpha / 2;
                        lmb = 1;
                        continue
                    end

                end

                if bad_iterations > 20
                    relative = 0;

                else
                end

                if (((relative < j_rm)) && all((sum(abs(new_defects), 2)) < obj.defects_max))
                    if obj.plot_end
                        obj.plot_xu([new_state_array; defects], new_u, time_array, obj.names, state_d, obj.order, "final")
                    end
                    return
                end

            end

        end %function



        %% Cost Function
        function J = cost(obj,state_array, state_d, u)
            [~, obj.horizon_disc] = size(state_array);
            J = 0;

            for i = 1:obj.horizon_disc - 1
                J = J + (state_array(:, i) - state_d(:, i))' * obj.Q * (state_array(:, i) - state_d(:, i)) + u(:,i)' * obj.R * u(:,i);
            end

            i = obj.horizon_disc;
            J = J + (state_array(:, i) - state_d(:, i))' * obj.Qn * (state_array(:, i) - state_d(:, i));
        end %function

        %% Backward
        function [L, l, A_, B_] = backward(obj,n_states, defects, x, xd, u)
            [n_controls, ~] = size(u);
            S = zeros(n_states, n_states, obj.horizon_disc);
            s = zeros(n_states, obj.horizon_disc); %deep horizon+1 and hight is n_states
            L = zeros(n_controls, n_states, obj.horizon_disc-1); %size depends both from the number of controls and states
            l = zeros(n_controls, obj.horizon_disc-1); % size depends from the number of controls
            %Fill the S and s matrix
            s(:, obj.horizon_disc) = obj.Qn * (x(:, obj.horizon_disc) - xd(:, obj.horizon_disc));
            S(:, :, obj.horizon_disc) = obj.Qn;
            A_ = zeros(n_states, n_states, (obj.horizon_disc - 1));
            B_ = zeros(n_states, n_controls, obj.horizon_disc - 1);

            %back
            for n = (obj.horizon_disc - 1):-1:1
                %N = (n_states * n - (n_states-1)):(n_states * n);
                [A_(:, :, n), B_(:, :, n)] = obj.model.linearization_discretization(u(:, n), x(:, n), 1);
                A = A_(:, :, n);
                B = B_(:, :, n);
                %compute r,h,G,H to simplify S and s computations
                P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
                r = obj.R * u(:, n); % the derivative of uRu
                q = obj.Q * (x(:, n) - xd(:, n));
                h = r + B.' * (s(:, n + 1) + S(:, :, n + 1) * (defects(:, n))); %gu + S(:, N + 4) * defects(:, n+1)
                G = P + B.' * S(:, :, n + 1) * A; %Gux
                H = obj.R + B.' * S(:, :, n + 1) * B; %Guu
                L(:, :, n) = -pinv(H) * G; %K
                l(:, n) = -pinv(H) * h; %d

                Gxx = obj.Q + A.' * S(:, :, n + 1) * A;
                gx = q + A' * (s(:, n + 1) + S(:, :, n + 1) * (defects(:, n)));

                % Gxu = G';
                % Gux = G;

                S(:, :, n) = Gxx - L(:, :, n).' * H * L(:, :, n);
                s(:, n) = gx + G.' * l(:, n) + L(:, :, n).' * (h + H * l(:, n));

            end

        end %function

        %% Forward Shoot
        function [x, u] = forward_shoot(obj,ix,  u, L, l, x_old)
            n_states = length(ix);
            x = repmat(ix, 1, obj.horizon_disc);
            for n = 1:obj.horizon_disc - 1
                u(:, n) = u(:, n) + l(:, n) + L(:, :,n) * (x(:, n) - x_old(:, n));
                x(:, n + 1) = obj.dynamics_rk4(x(:, n), u(:,n));
            end
        end %function

        %% Forward Linear Shoot
        function [x, u] = forward_linear_shoot(obj, x_old, u, L, l, defects, A_, B_, ~)
            x = x_old;
            [n_states, ~] = size(x); %#ok

            %x(:, 1) = ix; %initial state

            for n = 1:obj.horizon_disc - 1
                %N = (n_states * n - (n_states - 1)):(n_states * n);
                traj = 1; %
                contr= 0;%all(defects(:, n)==0);
                u(:, n) = u(:, n) + contr*(l(:, n) + L(:, :,n) * (x(:, n) - x_old(:, n)));
                %                 if J < 1e8
                A = A_(:, :,n);
                B = B_(:, :,n);
                x(:, n + 1) = x_old(:, n + 1) + (defects(:, n)) ... % take the previous value
                    +traj * ((A + B * L(:, :,n)) * (x(:, n) - x_old(:, n))+B * l(:, n));

                %                 else
                %                     x(:, n + 1) = obj.dynamics_euler(x(:, n), u(:, n), obj.dt);
                %                 end
            end

        end %function

        %% Forward Multi Shoot

        function [x, defects, u] = forward_multi_shoot(obj, x, u, x_old, L, l)

            %closedloop=0;

            n_states = length(x(:, 1));
            n_controls = length(u(:, 1));
            if obj.horizon_disc < 3
                obj.pieces = 1;
                k = 1;
            else
                %k=2.^(0:-0.1:(-pieces+1)/10);
                k = ones(1, obj.pieces) * 1;
            end

            [~, len] = size(x);
            len = floor(len / obj.pieces);

            x_start = zeros(n_states, obj.pieces);
            x_arrive = x_start;

            for i = 1:obj.pieces
                x_start(:, i) = x(:, len * (i - 1) + 1);
            end

            statess = zeros(n_states, obj.horizon_disc, obj.pieces);
            uss = zeros(n_controls, obj.horizon_disc - 1, obj.pieces);

            parfor i = 1:obj.pieces

                stato = zeros(n_states, obj.horizon_disc);
                us = zeros(n_controls, obj.horizon_disc - 1);

                t = len * (i - 1) + 1;
                t = t + (t == 0); %skip if t==0
                inizio = t;

                stato(:, inizio) = x_start(:, i);

                if i == obj.pieces %if is last piece do till the end
                    fine = obj.horizon_disc - 1;
                else
                    fine = len * (i);
                end
                %len(i-1)+1:len(i):
                for elem = inizio:fine
                    n = t;
                    %N = (n_states * n - (n_states - 1)):(n_states * n);
                    us(:, n) = u(:, n) + (k(i) * l(:, n) + L(:,:,n) * (stato(:, t) - x_old(:, n)));
                    stato(:, t + 1) = obj.dynamics_rk4(stato(:, t), us(:, elem));
                    t = t + 1;
                end

                x_arrive(:, i) = stato(:, t);

                if fine ~= (obj.horizon_disc - 1)
                    stato(:, t) = stato(:, t) * 0;
                end

                uss(:, :, i) = us(:, :);
                statess(:, :, i) = stato(:, :);
            end

            x = zeros(n_states, obj.horizon_disc);
            u = zeros(n_controls, obj.horizon_disc - 1);

            for i = 1:obj.pieces
                x(:, :) = x(:, :) + statess(:, :, i);
                u(:, :) = u(:, :) + uss(:, :, i);
            end

            defects = zeros(n_states, obj.horizon_disc);

            for i = 1:1:obj.pieces - 1
                defects(:, len*i+1) = (x_arrive(:, i) - x_start(:, i + 1));
            end
            %defects(:, end) = x_arrive(:, end) - state_d(:, end);
            defects = circshift(defects, -1, 2);
        end %function

        %% Plot
    end
    methods(Static)
        function plot_xu(x, u, time_array, names, xd, order, type)

            %plot_xu(x, u, time_array, names,xd,order)
            if type == ""
                return
            end

            if nargin < 3
                error("missing entries")
            end

            [h, ~] = size(x);

            if nargin < 4 || isempty(names) || length(names)<h
                init=length(names)+1;
                names =[names,init:h]; %repmat("", 1, h);
            end


            desired = true;

            if nargin < 5
                desired = false;
            end

            if nargin < 6||isempty(order)

                order = (1:h)';
            end

            if nargin < 7
                type = 'figure';
            end



            [h, l] = size(order);

            tiledlayout(h + 1, 1);

            for i = 1:h
                nexttile
                j = 1;
                lgd = [];

                while j <= l && ~isnan(order(i, j))

                    idx = order(i, j);
                    plot(time_array, x(idx, :))
                    hold on
                    lgd = [lgd, names(idx)+""]; %#ok

                    if desired && length(xd(:, 1)) >= idx
                        plot(time_array, xd(idx, :), LineStyle = "--")
                        hold on
                        lgd = [lgd, names(idx) + "_d"]; %#ok
                    end

                    title(type)

                    j =j+1;




                end
                legend(lgd,'Location','northeastoutside'); %interpreter =latex
                ylim('padded')
                xlim([-inf time_array(end)]);
                grid minor


            end
            xl = xlim;

            nexttile
            plot(time_array(1:end - 1), u)
            %legend("controls",'Location','northeastoutside')
            title(type+"   controls")
            xlim(xl);
            grid minor


            pause(0)
            set(gcf, 'Name', type, 'NumberTitle', 'off')
        end %function


    end

end