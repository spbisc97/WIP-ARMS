classdef (Abstract) Model
    methods
        function obj=Model()
            
        end
    end
    methods (Abstract , Static)
        dy=ForwardDynamics(y,u)
        [A, B] = linearization_discretization(u, y, discrete)
        new_state = euler_integration_fun(state, dy, dt)
    end
end