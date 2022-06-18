clear all, close all, 

syms x theta Dtheta Dx
syms DDx DDtheta  F mp mc l g 

g=10;  mp=1; l=2;  mc=5;


s = 1; % pendulum up (s=1)
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
disp(A)
disp(B)
eig(A)
pause
A=double(A);
B=double(B);

Q = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
R = .0001;

%%
det(ctrb(A,B))

%%
K = lqr(A,B,Q,R);

disp("start Int")
tspan = 0:.001:10;
if(s==-1)
    y0 = [0; 0; 0.1; 0];
    [t,y] = ode45(@(t,y)removecartpend(y,mp,mc,l,g,0,-K*(y-[0; 0; 0; 0])),tspan,y0);
elseif(s==1)
    y0 = [-3; 0; 0; 0];
% % [t,y] = ode45(@(t,y)((A-B*K)*(y-[0; 0; pi; 0])),tspan,y0);
    [t,y] = ode45(@(t,y)removecartpend(y,mp,mc,l,g,0,-K*(y-[0; 0; 0; 0])),tspan,y0);
else
    
end

for k=1:100:length(t)
    removedrawcartpend(y(k,:),mp,mc,l);
end

% function dy = pendcart(y,m,M,L,g,d,u)