function [A, B] = linearization_discretization(u, state, discrete)
    [A, B]=Twip.linearization_discretization(u,[0; 0; 0; 0; 0; 0],discrete);
end
