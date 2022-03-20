
clc;
clear all;
syms x theta Dtheta Dx
syms DDx DDtheta  F mp mc l g 

Q=diag([1,1,10,100]);
R=.001;

g=9.81;  mp=0.5; l=1;  mc=3;

eq(1)= -mp*g*sin(theta)==mp*DDx*cos(theta)-mp*l*DDtheta;
eq(2)= F+mp*l*DDtheta*cos(theta)-mp*l*Dtheta^2*sin(theta)==(mp+mc)*DDx;
sol=solve(eq,[DDx,DDtheta]);



A=sym(zeros(4,4));
A(1,2)=1;
A(3,4)=1;
A(2,:)=subs(jacobian(sol.DDx,[x Dx theta Dtheta]),[theta Dtheta],[0,0]);
A(4,:)=subs(jacobian(sol.DDtheta,[x Dx theta Dtheta]),[theta Dtheta],[0,0]);

B=sym(zeros(4,1));
B(2)=subs(jacobian(sol.DDx,F),[theta Dtheta],[0,0]);
B(4)=subs(jacobian(sol.DDtheta,F),[theta Dtheta],[0,0]);

% lqr(A,B,Q,R)
