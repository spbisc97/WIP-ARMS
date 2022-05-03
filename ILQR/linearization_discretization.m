function [A,B]=linearization_discretization(u,state)
Ts=0.01;
A=1+2*state;
B=1;

B = pinv(eye(1) - A*Ts/2.)*B*sqrt(Ts);
A = (eye(1) + A*Ts/2.)*pinv(eye(1) - A*Ts/2.);


end
