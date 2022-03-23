



clear all, close all, clc

tspan = 0:.0005:10;

y0 = [0.9; 0; 0; 0;0;0];
[t,y] = ode45(@(t,y)3dwip(y,mp,mc,l,g,0,r,0),tspan,y0);



syms x theta Dtheta Dx
syms DDx DDtheta  F mp mc l g 


syms DDpsi DDtheta theta Dtheta psi
syms Tr r mw g mp l Te  tau


% g=10;  l=2;  mc=5;
% mp=2; mw=1; r=1;
% tau=F;










A=sym(zeros(12,12));
A(1,2)=1;
A(2,:)=subs(jacobian(sol.DDpsi_,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);
A(3,4)=1;
A(4,:)=subs(jacobian(sol.DDpsi,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);
A(5,6)=1;
A(6,:)=subs(jacobian(sol.DDtheta,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);
% A(7,8)=1
% A(8,:)=subs(jacobian(sol.DDtheta,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);
% A(9,10)=1
% A(10,:)=subs(jacobian(sol.DDtheta,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);
% A(11,12)=1
% A(12,:)=subs(jacobian(sol.DDtheta,[x Dx theta Dtheta]),[theta Dtheta],[stab,0]);

B=sym(zeros(6,1));
B(1)=subs(jacobian(sol.DDpsi,F),[theta Dtheta],[stab,0]);
B(2)=
B(3)=

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
R = .001;

%%
det(ctrb(A,B))

%%
K = [0,0,0,0]%lqr(A,B,Q,R);
s=-1;
disp("start Int")
tspan = 0:.0005:10;
if(s==-1)
    y0 = [0; 0; stab+0.9; 0];
    [t,y] = ode45(@(t,y)3dwip(y,mp,mc,l,g,0,r,K*(y-[0; 0; stab; 0])),tspan,y0);
elseif(s==1)
    y0 = [-3; 0; 0; 0];
% % [t,y] = ode45(@(t,y)((A-B*K)*(y-[0; 0; pi; 0])),tspan,y0);
    [t,y] = ode45(@(t,y)3dwip(y,mp,mc,l,g,0,r,-K*(y-[0; 0; 0; 0])),tspan,y0);
else
    
end

for k=1:100:length(t)
    drawwip(y(k,:),mp,mc,l,r);
end

%% 

function drawwip(y,m,M,L,r)
x = y(1);
th = -y(3)+pi;

% kinematics
% x = 3;        % cart position
% th = 3*pi/2;   % pendulum angle

% dimensions
% L = 2;  % pendulum length
W = r;  % cart width
H = r; % cart height
wr = 0; % wheel radius
mr = .3*sqrt(m); % mass radius

% positions
% y = wr/2; % cart vertical position
y = wr/2+H/2; % cart vertical position


px = x + L*sin(th);
py = y - L*cos(th);

plot([-10 10],[0 0],'w','LineWidth',2)
hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1])


plot([x px],[y py],'w','LineWidth',2)

rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1])

% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim([-5 5]);
ylim([-2 2.5]);
set(gca,'Color','k','XColor','w','YColor','w')
set(gcf,'Position',[10 900 800 400])
set(gcf,'Color','k')
set(gcf,'InvertHardcopy','off')   

% box off
drawnow
hold off

end
%% 

function dy = d3wip(y,mp,mc,l,g,d,r,u_r,u_l)
theta=y(3);
dtheta=y(4);
psi_l
dpsi_l
psi_r
dpsi_r



dy(1,1) = y(2);
dy(2,1) = (5.0000e-04*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40);
dy(3,1) = y(4);
dy(4,1) = -(0.0080*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);
dy(5,1) = y(6);
dy(6,1) = -(0.0080*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);

disp(dy.')
%pause
end