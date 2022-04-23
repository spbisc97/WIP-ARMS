function [A,B]=linearization_discretization(u,state)
A=1+2*state;
B=1;

B = pinv(eye(1) - A*0.01/2.)*B*sqrt(0.01);
A = (eye(1) + A*0.01/2.)*pinv(eye(1) - A*0.01/2.);
end
