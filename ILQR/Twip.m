classdef Twip < Model

    methods (Static)

        function dy = ForwardDynamics(y, u)

            %ForwardDynamics(~,y,u)

            dy = zeros(6, 1);

            T_l = u(1);
            T_r = u(2);

            M = 1;
            l = 2;
            g = 10;
            r = 1;
            d = 2;
            I_W = 0.5;
            M_W = 0.3;
            I_B = 0.5;
            I_W_d = 0.3;
            I_Z = 0.3;

            phi = y(1);
            phi_dot = y(2);
            x = y(3);
            x_dot = y(4);
            theta = y(5);
            theta_dot = y(6);

            u_2 = T_l + T_r;
            u_1 = T_l - T_r;

            Delta = -M^2 * l^2 * r^2 -2 * M * l^2 * I_W - 2 * M * l^2 * M_W * r^2 - I_B * M * r^2 - 2 * I_B * I_W - 2 * M_W * r^2 * I_B;
            Omega = -M_W * d^2 * r^2 - 4 * I_W_d * r^2 - 2 * I_Z * r^2 - d^2 * I_W - 2 * M * l^2 * r^2;

            num = -M * l * (theta_dot^2 * l * cos(phi) + g) * (M * r^2 * sin(phi) + 2 * I_W * sin(phi) + 2 * M_W * r^2 * sin(phi)) + M^2 * l^2 * r^2 * phi_dot^2 * sin(phi) * cos(phi);
            den = Delta + M^2 * l^2 * r^2 * cos(phi)^2;
            dy(2) = num / den + ((M * r^2 + 2 * I_W + 2 * M_W * r^2 + M * l * r * cos(phi)) * u_2) / (Delta + M^2 * l^2 * r^2 * cos(phi)^2);

            num = M * l * r^2 * (M * l * sin(phi) * cos(phi) * (theta_dot^2 * l * cos(phi) + g) - phi_dot^2 * sin(phi) * M * l^2 - phi_dot^2 * sin(phi) * I_B) - (M * l * r^2 * cos(phi) + M * r * l^2 + I_B * r) * u_2;
            den = Delta + M^2 * l^2 * r^2 * cos(phi)^2;
            dy(4) = num / den;

            num = 2 * M * r^2 * l^2 * theta_dot * phi_dot * sin(2 * phi) + r * d * u_1;
            den = Omega + 2 * M * l^2 * r^2 * cos(phi)^2;
            dy(6) = num / den;

            dy(1) = phi_dot;

            dy(3) = x_dot;

            dy(5) = theta_dot;

        end

        function [A_jac, B_jac] = auxiliary_jacobian(symbolic)

            if nargin < 1
                symbolic = false;
            end

            syms theta theta_dot theta_ddot
            syms phi phi_dot phi_ddot
            syms x x_dot x_ddot

            syms T_l;
            syms T_r;

            if symbolic
                syms M;
                syms l;
                syms g;
                syms r;
                syms d;
                syms I_W;
                syms M_W;
                syms I_B;
                syms I_W_d;
                syms I_Z;

            else
                M = 1;
                l = 2;
                g = 10;
                r = 1;
                d = 2;
                I_W = 0.5;
                M_W = 0.3;
                I_B = 0.5;
                I_W_d = 0.3;
                I_Z = 0.3;
            end

            u_2 = T_l + T_r;
            u_1 = T_l - T_r;

            Delta = -M^2 * l^2 * r^2 -2 * M * l^2 * I_W - 2 * M * l^2 * M_W * r^2 - I_B * M * r^2 - 2 * I_B * I_W - 2 * M_W * r^2 * I_B;
            Omega = -M_W * d^2 * r^2 - 4 * I_W_d * r^2 - 2 * I_Z * r^2 - d^2 * I_W - 2 * M * l^2 * r^2;

            num = -M * l * (theta_dot^2 * l * cos(phi) + g) * (M * r^2 * sin(phi) + 2 * I_W * sin(phi) + 2 * M_W * r^2 * sin(phi)) + M^2 * l^2 * r^2 * phi_dot^2 * sin(phi) * cos(phi);
            den = Delta + M^2 * l^2 * r^2 * cos(phi)^2;
            dy(2) = num / den + ((M * r^2 + 2 * I_W + 2 * M_W * r^2 + M * l * r * cos(phi)) * u_2) / (Delta + M^2 * l^2 * r^2 * cos(phi)^2);

            num = M * l * r^2 * (M * l * sin(phi) * cos(phi) * (theta_dot^2 * l * cos(phi) + g) - phi_dot^2 * sin(phi) * M * l^2 - phi_dot^2 * sin(phi) * I_B) - (M * l * r^2 * cos(phi) + M * r * l^2 + I_B * r) * u_2;
            den = Delta + M^2 * l^2 * r^2 * cos(phi)^2;
            dy(4) = num / den;

            num = 2 * M * r^2 * l^2 * theta_dot * phi_dot * sin(2 * phi) + r * d * u_1;
            den = Omega + 2 * M * l^2 * r^2 * cos(phi)^2;
            dy(6) = num / den;

            dy(1) = phi_dot;

            dy(3) = x_dot;

            dy(5) = theta_dot;

            dy = dy.';

            A_jac = jacobian(dy, [phi phi_dot x x_dot theta theta_dot]);
            B_jac = jacobian(dy, [T_l, T_r]);

        end

        function [A, B] = linearization_discretization(u, y, discrete, A_jac, B_jac)

            y=y(:);
            u = u(:);

            if nargin < 4

                phi = y(1);
                phi_dot = y(2);
                x = y(3);
                x_dot = y(4);
                theta = y(5);
                theta_dot = y(6);
                T_l = u(1);
                T_r = u(2);

                A = [[0, 1, 0, 0, 0, 0];
                    [(4 * phi_dot^2 * cos(phi)^2 - 4 * phi_dot^2 * sin(phi)^2 + (52 * theta_dot^2 * sin(phi)^2) / 5 - (13 * cos(phi) * (4 * cos(phi) * theta_dot^2 + 20)) / 5) / (4 * cos(phi)^2 - 117/10) - (2 * sin(phi) * (T_l + T_r)) / (4 * cos(phi)^2 - 117/10) - (8 * cos(phi) * sin(phi) * (- 4 * cos(phi) * sin(phi) * phi_dot^2 + (13 * sin(phi) * (4 * cos(phi) * theta_dot^2 + 20)) / 5)) / (4 * cos(phi)^2 - 117/10)^2 + (8 * cos(phi) * sin(phi) * (T_l + T_r) * (2 * cos(phi) + 13/5)) / (4 * cos(phi)^2 - 117/10)^2, (8 * phi_dot * cos(phi) * sin(phi)) / (4 * cos(phi)^2 - 117/10), 0, 0, 0, - (104 * theta_dot * cos(phi) * sin(phi)) / (5 * (4 * cos(phi)^2 - 117/10))];
                    [0, 0, 0, 1, 0, 0];
                    [- (9 * phi_dot^2 * cos(phi) + 4 * sin(phi)^2 * (2 * cos(phi) * theta_dot^2 + 10) - 2 * sin(phi) * (T_l + T_r) - 4 * cos(phi)^2 * (2 * cos(phi) * theta_dot^2 + 10) + 8 * theta_dot^2 * cos(phi) * sin(phi)^2) / (4 * cos(phi)^2 - 117/10) - (8 * cos(phi) * sin(phi) * (9 * sin(phi) * phi_dot^2 + (T_l + T_r) * (2 * cos(phi) + 9/2) - 4 * cos(phi) * sin(phi) * (2 * cos(phi) * theta_dot^2 + 10))) / (4 * cos(phi)^2 - 117/10)^2, - (18 * phi_dot * sin(phi)) / (4 * cos(phi)^2 - 117/10), 0, 0, 0, (16 * theta_dot * cos(phi)^2 * sin(phi)) / (4 * cos(phi)^2 - 117/10)];
                    [0, 0, 0, 0, 0, 1];
                    [(16 * cos(phi) * sin(phi) * (2 * T_l - 2 * T_r + 8 * phi_dot * theta_dot * sin(2 * phi))) / (8 * cos(phi)^2 - 13)^2 + (16 * phi_dot * theta_dot * cos(2 * phi)) / (8 * cos(phi)^2 - 13), (8 * theta_dot * sin(2 * phi)) / (8 * cos(phi)^2 - 13), 0, 0, 0, (8 * phi_dot * sin(2 * phi)) / (8 * cos(phi)^2 - 13)]; ];

                B = [[0, 0];
                    [(2 * cos(phi) + 13/5) / (4 * cos(phi)^2 - 117/10), (2 * cos(phi) + 13/5) / (4 * cos(phi)^2 - 117/10)];
                    [0, 0];
                    [- (2 * cos(phi) + 9/2) / (4 * cos(phi)^2 - 117/10), - (2 * cos(phi) + 9/2) / (4 * cos(phi)^2 - 117/10)];
                    [0, 0];
                    [2 / (8 * cos(phi)^2 - 13), -2 / (8 * cos(phi)^2 - 13)]; ];

            else
                syms theta theta_dot
                syms phi phi_dot
                syms x x_dot
                syms T_l;
                syms T_r;
                vars = [phi; phi_dot; x; x_dot; theta; theta_dot; T_l; T_r];
                values = [y; u];
                A = subs(A_jac, vars, values);
                B = subs(B_jac, vars, values);
            end

            if discrete == 1
                Ts = 0.01;
                n = 6;
                Ad = (eye(n) + A * Ts / 2.) * pinv(eye(n) - A * Ts / 2.); %tustin, bilinear trans
                B = pinv(eye(n) - A * Ts / 2.) * B * sqrt(Ts);
                Bd = pinv(A) * (Ad - eye(n)) * B;
                A = Ad;
                B = Bd;
            end

            A = double(A);
            B = double(B);
        end

        function new_state = euler_integration_fun(state, dy, dt)
            %integrate as euler
            %%test with single value
            %new_state=state+dy*dt;
            %%test for inverted pendulum
            new_state(2) = state(2) + dy(2) * dt;
            new_state(1) = state(1) + new_state(2) * dt;

            new_state(4) = state(4) + dy(4) * dt;
            new_state(3) = state(3) + new_state(4) * dt;

            new_state(6) = state(6) + dy(6) * dt;
            new_state(5) = state(5) + new_state(5) * dt;
            new_state = new_state.';
        end
    end

end
