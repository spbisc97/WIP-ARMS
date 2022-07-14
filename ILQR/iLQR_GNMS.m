classdef iLQR_GNMS
    properties
        %running cost multiplier
        Q  = 1
        %control cost multiplier
        R = 0.1
        %final cost multiplier
        Qn  =1
        %number of pieces for multishooting, equal to one when
        %using SS functions
        pieces  =16
        %defects max, can be set as single double or a vertical
        %array with same lenght of the states
        defects_max =0.001
        %prevision horizon
        horizon double = 10;
        %horzon discretized, the steps number, automatically
        %updated when defects max is changed
        horizon_disc=1001
        %the order for printing states, same line same tile, fill
        %empty with nan: [1,2,nan,nan;3,4,5,6]
        order=[]
        %the names of the variables ordered as x
        names=[]
        %system model to optimize
        model
        %how often print the iteration,1 for each time
        plot_steps  = 1;
        %plot the initial start without control
        plot_start  =true ;
        %plot the end configuration
        plot_end = true ;
        %dt of optimizer, needs to be the same of the model
        dt=0.01
        %pause after each plot, inf to need button press to skip
        plot_duration=0;
        %optional: figure where the plots will be drown
        plot_figure
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
                obj.Qn=Q;
            end

            if nargin > 3
                obj.Qn=Qn;
            end
            if nargin > 4
                obj.pieces=pieces;
            end
            if nargin > 5
                obj.defects_max=defects_max;
            end
        end %function

        %%main functions
        function u = MS_iLQR(obj,istate, state_d, it, u)
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
            u = [u(:, 2:end), repmat(u(:,end),1, (obj.horizon_disc - 1) - (usz - 1))];
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
                obj.plot_xu(state_array, u, time_array, obj.names,state_d,defects,obj.order,"multi start",obj.plot_duration,obj.plot_figure)
            end
            J = obj.cost(state_array, state_d, new_u) * 1e3;
            new_J = J;
            disp("J")
            disp(J)

            %start the optimizing iterations
            for iteration = 1:iterations - 1
                [L, l, A_, B_] = backward(obj,n_states,  defects, state_array, state_d, u);
                
                [new_state_array, new_u] = forward_linear_shoot(obj, state_array,state_array, u,lmb * L, alpha * l, defects, A_, B_, new_J);

                if mod(iteration, obj.plot_steps) == 0

                    obj.plot_xu( new_state_array, new_u, time_array, obj.names, state_d,defects, obj.order, "linear",obj.plot_duration,obj.plot_figure)
                end


                [new_state_array, new_defects, new_u] = obj.forward_multi_shoot(new_state_array, new_u, state_array, lmb*L, alpha*l);

                %plot before exit to understand what is happening
                if mod(iteration, obj.plot_steps) == 0
                    obj.plot_xu(new_state_array, new_u, time_array,obj.names, state_d,new_defects,obj.order, "multi",obj.plot_duration,obj.plot_figure)

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
                        obj.plot_xu(new_state_array, new_u, time_array, obj.names, state_d, defects,obj.order, "final",obj.plot_duration,obj.plot_figure)
                    end
                    return
                end

            end

        end %function

        function u = SS_iLQR(obj,istate, state_d, it, u)
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
            u = [u(:, 2:end), repmat(u(:,end),1, (obj.horizon_disc - 1) - (usz - 1))];
            % u=u*0; reset suggested controls
            state_array = ones(1, obj.horizon_disc) .* istate;
            new_u = u;
            alpha = 1;
            lmb = 1;
            time_array = it:obj.dt:(obj.horizon + it);

            % !finished initialization


            L = zeros(n_controls,n_states, obj.horizon_disc ); %size depends both from the number of controls and states
            l = zeros(n_controls, obj.horizon_disc); % size depends from the number of controls
            [state_array, u] = forward_shoot(obj, state_array, u, state_array, L, l);
            if obj.plot_start
                obj.plot_xu(state_array, u, time_array, obj.names,state_d,[],obj.order,"start",obj.plot_duration,obj.plot_figure)
            end
            J = obj.cost(state_array, state_d, new_u);
            new_J = J;
            new_state_array=state_array;
            disp("J")
            disp(J)

            %start the optimizing iterations
            for iteration = 1:iterations - 1
                [L, l, ~, ~] = backward_ilqr_ss(obj,n_states, state_array, state_d, u);

                [new_state_array, new_u] = obj.forward_shoot(new_state_array, new_u, state_array, lmb*L, alpha*l);

                %plot before exit to understand what is happening
                if mod(iteration, obj.plot_steps) == 0
                    obj.plot_xu(new_state_array, new_u, time_array,obj.names, state_d,[],obj.order, "ss ilqr",obj.plot_duration,obj.plot_figure)

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

                T = table(new_J, J, relative, bad_iterations, alpha, 'VariableNames', {'new_cost', 'prev_cost', 'relative', 'last bad', 'alpha'});

                disp(T)

                if bad_iterations == 0
                    state_array = new_state_array;
                    J = new_J;
                    u = new_u;
                    alpha = 1;
                    lmb = 1;
                else

                    if ~isnan(new_J) && ~isinf(new_J)
                        alpha = alpha / 2;
                        lmb = 1;
                        state_array = new_state_array;
                        J = new_J;
                        u = new_u;
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

                if (((relative < j_rm)))
                    if obj.plot_end
                        obj.plot_xu(new_state_array, new_u, time_array, obj.names, state_d, [],obj.order, "final",obj.plot_duration,obj.plot_figure)
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
            if isempty(defects)
                defects=zeros(n_states, obj.horizon_disc);
            end
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
                %check 
                %disp(A*x(:,n)+B*u(:,n) -x(:,n+1))%error on discretization
                %compute r,h,G,H to simplify S and s computations
                P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
                
                %R and Q derivatives
                r = obj.R * u(:, n); % the derivative of uRu
                q = obj.Q * (x(:, n) - xd(:, n));

                %gu calculation
                h = r + B.' * (s(:, n + 1) + S(:, :, n + 1) * (defects(:, n)));
                %Gux
                G = P + B.' * S(:, :, n + 1) * A; 
                %Guu
                H = obj.R + B.' * S(:, :, n + 1) * B;
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

        function [L, l, A_, B_] = backward_ilqr_ss(obj,n_states, x, xd, u)
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
                %check 
                %compute r,h,G,H to simplify S and s computations
                P = 0; %repmat([0], [1, n_states]); %set mixed weight to zero delta_u*P*delta_x
                
                %R and Q derivatives
                r = obj.R * u(:, n); % the derivative of uRu
                q = obj.Q * (x(:, n) - xd(:, n));

                %gu calculation
                gu = r + B.' * (s(:, n + 1));

                gx = q + A' * (s(:, n + 1));

                %Gux
                Gxu =P+ A.' * S(:, :,n + 1) * B;
                Gux =P+ B.' * S(:, :,n + 1) * A;
                %Guu
                Guu = obj.R + B.' * S(:, :, n + 1) * B;
                L(:, :, n) = -pinv(Guu) * Gux; %K
                l(:, n) = -pinv(Guu) * gu; %d

                Gxx = obj.Q + A.' * S(:, :, n + 1) * A;


                S(:,:, n) = Gxx + L(:, :,n).' * Guu * L(:,:,n) + Gxu * L(:, :,n) + L(:, :,n).' * Gux;
                s(:, n) = gx + L(:, :,n).' * gu + L(:, :,n).' * Guu * l(:, n) + Gxu * l(:, n);
    
            end

        end %function

        %% Forward Shoot
        function [x, u] = forward_shoot(obj,x,  u, x_old,L, l)
            for n = 1:obj.horizon_disc - 1
                u(:, n) = u(:, n) + l(:, n) + L(:, :,n) * (x(:, n) - x_old(:, n));
                x(:, n + 1) = obj.dynamics_rk4(x(:, n), u(:,n));
            end
        end %function

        %% Forward Linear Shoot
        function [x, u] = forward_linear_shoot(obj,x, x_old, u, L, l, defects, A_, B_, ~)
            [n_states, ~] = size(x); %#ok
            for n = 1:obj.horizon_disc - 1
                contr=0;
                traj=1;
                u(:, n) = u(:, n) + contr*(l(:, n) + L(:, :,n) * (x(:, n) - x_old(:, n)));
                A = A_(:, :,n);
                B = B_(:, :,n);
                x(:, n + 1) = x_old(:, n + 1) + (defects(:, n)) ... % starting point
                    +traj * ((A + B * L(:, :,n)) * (x(:, n) - x_old(:, n))+B * l(:, n));%adjustment


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
            if coder.target('MATLAB')
                for i = 1:obj.pieces
                    [ x_arrive(:, i) ,uss(:, :, i) ,statess(:, :, i)]=multi_aux(obj,i,x_start,len,x_old,u,l,L,n_states,n_controls);
                end
            else
                parfor i = 1:obj.pieces
                [ x_arrive(:, i) ,uss(:, :, i) ,statess(:, :, i)]=multi_aux(obj,i,x_start,len,x_old,u,l,L,n_states,n_controls);
                end
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


        function [x_arrive,us,stato]=multi_aux(obj,i,x_start,len,x_old,u,l,L,n_states,n_controls)
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
                %N = (n_states * n - (n_states - 1)):(n_states * n);
                us(:, t) = u(:, t) + (l(:, t) + L(:,:,t) * (stato(:, t) - x_old(:, t)));
                stato(:, t + 1) = obj.dynamics_rk4(stato(:, t), us(:, t));
                t = t + 1;
            end

            x_arrive= stato(:, t);

            if fine ~= (obj.horizon_disc - 1)
                stato(:, t) = stato(:, t) * 0;
            end

            
        end
    end
    methods(Static)
        function plot_xu(x, u, time_array, names, xd,defects,order, type,plot_duration,plot_figure)

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

            if nargin < 7||isempty(order)

                order = (1:h)';
            end

            if nargin < 8
                type = 'GNMS_iLQR';
            end
            if nargin < 9
                plot_duration = inf;
            end
            if nargin < 10 || isempty(plot_figure)% || class(plot_figure)~="matlab.ui.Figure"
                figure_name="GNMS_iLQR";
            else
                set(0,'currentfigure',plot_figure)
                figure_name=plot_figure.Name;
            end

            
            [h, l] = size(order);

            if isempty(defects)
                tiles_n=h+1;
                def=false;
            else
                tiles_n=h+2;
                def=true;
                defplot  = double(any(circshift(defects,1,2)~=0,1));
                defplot(defplot==0)=nan;
            end





            tiledlayout(tiles_n, 1,'Padding','Compact','Tilespacing','Compact');

            for i = 1:h
                nexttile(i)
                j = 1;
                lgd = [];

                while j <= l && ~isnan(order(i, j))

                    idx = order(i, j);
                    plot(time_array, x(idx, :))
                    hold on
                    lgd = [lgd, names(idx)+""]; %#ok

                    if def
                        plot(time_array, x(idx, :).*defplot,"--o")
                        lgd = [lgd, ""]; %#ok
                    end
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
            lgd=[];
            if def
                nexttile(h+1)
                plot(time_array(1:end), defects)
                for i=names
                    lgd=[lgd,"def-"+i];
                end
                legend(lgd,'Location','northeastoutside')
                title(type+"  defects")
                ylim('padded')
                xlim(xl);
                grid minor
            end




            nexttile(tiles_n)
            plot(time_array(1:end - 1), u)
            %legend("controls",'Location','northeastoutside')
            title(type+"   controls")
            ylim('padded')
            xlim(xl);
            grid minor

            if isinf(plot_duration)
                pause
            else
                pause(plot_duration)
            end
            set(gcf, 'Name', figure_name, 'NumberTitle', 'off')
        end %function


    end

end
