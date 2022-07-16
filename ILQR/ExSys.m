classdef ExSys < Model
    properties
        dt=0.01
    end
    methods (Static)
        function new_state = euler_integration_fun(state, dy, dt)
            new_state(1) = state(1) + dy(1) * dt;
            new_state = new_state.';
        end
        function dy = ForwardDynamics(y, u)
            u=u(:);
            y=y(:);
            dy(1,1)=(1+y)*y+u; 
        end
        function [A, B] = linearization_discretization(~, state, discrete)
            if nargin < 3
                discrete = 0;
            end
            Ts=0.01;
            A=(1+2*state);
            B=1;
            if discrete == 1
                n = 1;
                Ad = (eye(n) + A * Ts / 2.) * pinv(eye(n) - A * Ts / 2.); %tustin, bilinear trans
                B = pinv(eye(n) - A * Ts / 2.) * B * sqrt(Ts);
                Bd = pinv(A) * (Ad - eye(n)) * B;
                A = Ad;
                B = Bd;
            end
        end
    end
end
