


clear all, close all, clc


syms x theta Dtheta Dx
syms DDx DDtheta  F mp mc l g 


syms DDpsi DDtheta theta Dtheta psi
syms Tr r mw g mp l Te  %tau


g=10;  mp=1; l=2;  mc=5;
mp=5; mw=1; r=1;
 tau=F;

stab=0;

eq(1)= tau/r+Te*sin(theta)==mw*DDpsi/r;
eq(1)= tau+Tr==r*mw*DDpsi; 
eq(2)= -Tr+sin(theta)*l*mp*g==mp*l^2/12*DDtheta;
eq(3)= -Te*sin(theta)==mp*(DDpsi/r+l*DDtheta*cos(theta)-l*Dtheta^2*sin(theta));
eq(4)= -Te*cos(theta)-mp*g==mp*(-l*DDtheta*sin(theta)-l*Dtheta^2*cos(theta));


sol=solve(eq,[DDtheta DDpsi Te Tr]);
sol.DDtheta=simplify(sol.DDtheta);
sol.DDpsi=simplify(sol.DDpsi);
disp(sol)


A=sym(zeros(4,4));
A(1,2)=1;
A(3,4)=1;
A(2,:)=subs(jacobian(sol.DDpsi,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);
A(4,:)=subs(jacobian(sol.DDtheta,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);

B=sym(zeros(4,1));
B(2)=subs(jacobian(sol.DDpsi,F),[theta Dtheta],[stab,0]);
B(4)=subs(jacobian(sol.DDtheta,F),[theta Dtheta],[stab,0]);
disp(A)
disp(B)
eig(A)
pause
A=double(A);
B=double(B);

Q = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 10];
R = .0001;

%%
det(ctrb(A,B))

%%
K =[0,0,0,0]% lqr(A,B,Q,R);
s=-1;
disp("start Int")
tspan = 0:.001:1;
if(s==-1)
    y0 = [0; 0; stab; 0];
    [t,y] = ode45(@(t,y)removewip(y,mp,mc,l,g,0,r,-K*(y-[0; 0; stab; 0])),tspan,y0);
elseif(s==1)
    y0 = [-3; 0; 0; 0];
% % [t,y] = ode45(@(t,y)((A-B*K)*(y-[0; 0; pi; 0])),tspan,y0);
    [t,y] = ode45(@(t,y)removewip(y,mp,mc,l,g,0,r,-K*(y-[0; 0; 0; 0])),tspan,y0);
else
    
end

for k=1:100:length(t)
    drawwip(y(k,:),mp,mc,l,r);
end

% function dy = pendcart(y,m,M,L,g,d,u)