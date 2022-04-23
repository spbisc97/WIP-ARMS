



 close all, clc





syms x theta Dtheta Dx
syms DDx DDtheta  F mp mc l g 
global WIP;
global Ag Bg;
WIP=WIP3dModel(9.81);


syms ddpsi_l ddtheta dtheta theta dpsi_l psi_l dpsi_r psi_r 
syms Tr r mw g mp l Te  tau u_l u_r

syms theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi
stab=0;

% g=10;  l=2;  mc=5;
% mp=2; mw=1; r=1;
% tau=F;


g=9.81;










A=sym(zeros(12,12));
A(1,2)=1;
A(2,:)=jacobian(WIP.ddtheta,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
A(3,4)=1;
A(4,:)=jacobian(WIP.ddpsi_l,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
A(5,6)=1;
A(6,:)=jacobian(WIP.ddpsi_r,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
A(7,8)=1;
A(8,:)=jacobian(WIP.ddx,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
A(9,10)=1;
A(10,:)=jacobian(WIP.ddy,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
A(11,12)=1;
A(12,:)=jacobian(WIP.ddphi,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
Ag=A;

A=subs(A,[theta dtheta psi_r dpsi_r psi_l dpsi_l u_l u_r],[stab,0,0,0,0,0,0,0]);
B=sym(zeros(6,2));
B(1)=0;
B(2,:)=(jacobian(WIP.ddtheta,[u_l,u_r]));
B(4,:)=(jacobian(WIP.ddpsi_l,[u_l,u_r]));
B(6,:)=(jacobian(WIP.ddpsi_l,[u_l u_r]));
B(8,:)=(jacobian(WIP.ddx,[u_l,u_r]));
B(10,:)=(jacobian(WIP.ddy,[u_l,u_r]));
B(12,:)=(jacobian(WIP.ddphi,[u_l u_r]));
Bg=B;
B=subs(B,[theta dtheta psi_r dpsi_r psi_l dpsi_l u_l u_r],[stab,0,0,0,0,0,0,0]);

disp(A)
disp(B)
eig(A)
pause
A=double(A);
B=double(B);

Q = diag([10,10,1,1,1,1]);
R = [0.01 0;
    0 0.01];



K=lqr_fun(A,B);
%disp(K)
pause
s=-1;
disp("start Int")
tspan = 0:.001:1;
if(s==-1)
    %y0 = [0.001; 0; 0; 0;0;0;0;0;0;0;0;0];
    %[t,y] = ode45(@(t,y)d3wip(y,0,0,0,9.81,0,0,-K*(y(1:2)-[0;0])),tspan,y0);
    %[t,y] = ode45(@(t,y)d3wip(y,0,0,0,9.81,0,0,-K*([y(1:2);y(7:10)]-[0;0;0;0;0;0])),tspan,y0);
    y0 = [0.001; 0; 2; 0;-2;0;0;0;0;0;12/10;0];
    [t,y] = ode45(@(t,y)d3wip(y,0,0,0,9.81,0,0,-lqr__fun(y)*([y(1:2);y(11:12)]-[0;0;0;0])),tspan,y0);

elseif(s==1)
    y0 = [2; 0; 0; 0;0;0;0;0;0;0;0;0];
% % [t,y] = ode45(@(t,y)((A-B*K)*(y-[0; 0; pi; 0])),tspan,y0);
    [t,y] = ode45(@(t,y)d3wip(y,mp,mc,l,g,0,r,[0;0]),tspan,y0);
else
    
end
figure('Name','body')
tiledlayout(1,2);
nexttile
plot(t,y(:,1).')
legend('theta')
nexttile
plot(t,y(:,2).')
legend('dtheta')

figure('Name','Wheels')
tiledlayout(2,2);
nexttile
plot(t,y(:,3).')
legend('psi_r')
nexttile
plot(t,y(:,4).')
legend('dpsi_r')
nexttile
plot(t,y(:,5).')
legend('psi_l')
nexttile
plot(t,y(:,6).')
legend('dpsi_l')

figure('Name','Se(2)')
tiledlayout(3,2);
nexttile
plot(t,y(:,7).')
legend('x')
nexttile
plot(t,y(:,8).')
legend('dx')
nexttile
plot(t,y(:,9).')
legend('y')
nexttile
plot(t,y(:,10).')
legend('dy')
nexttile
plot(t,y(:,11).')
legend('phi')
nexttile
plot(t,y(:,12).')
legend('dphi')




%legend('theta','dtheta','psi_r','dpsi_r','psi_l','dpsi_l','x','dx','y','dy','phi','dphi')

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

function dy = d3wip(y,~,~,~,g,~,~,u)
theta=y(1);
dtheta=y(2);
psi_r=y(3);
dpsi_r=y(4);
psi_l=y(5);
dpsi_l=y(6);
u_r=u(1);
u_l=u(2);




dy(1,1) = y(2);
dy(2,1) = (5.0000e-04*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40);
dy(3,1) = y(4);
dy(4,1) = -(0.0080*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);
dy(5,1) = y(6);
dy(6,1) = -(0.0080*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);
dy(7,1) = y(8);
dy(8,1) = (1.5000e-05*cos(0.3000*psi_l - 0.3000*psi_r)*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40) - (1.2000e-04*cos(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - (1.2000e-04*cos(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - 0.3000*dpsi_l*sin(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta) + 0.3000*dpsi_r*sin(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta);
dy(9,1) =y(10);
dy(10,1) = (1.2000e-04*sin(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) + (1.2000e-04*sin(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - (1.5000e-05*sin(0.3000*psi_l - 0.3000*psi_r)*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40) - 0.3000*dpsi_l*cos(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta) + 0.3000*dpsi_r*cos(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta);
dy(11,1) =y(11);
dy(12,1) =(0.0024*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - (0.0024*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);

end
function [A_,B_]=linWip(state)
%state=[stab,0,0,0,0,0,0,0]
global WIP Ag Bg
syms ddpsi_l ddtheta dtheta theta dpsi_l psi_l dpsi_r psi_r 
syms Tr r mw g mp l Te  tau u_l u_r

syms theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi

% A=sym(zeros(12,12));
% A(1,2)=1;
% A(2,:)=jacobian(WIP.ddtheta,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
% A(3,4)=1;
% A(4,:)=jacobian(WIP.ddpsi_l,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
% A(5,6)=1;
% A(6,:)=jacobian(WIP.ddpsi_r,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
% A(7,8)=1;
% A(8,:)=jacobian(WIP.ddx,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
% A(9,10)=1;
% A(10,:)=jacobian(WIP.ddy,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
% A(11,12)=1;
% A(12,:)=jacobian(WIP.ddphi,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi]);
% 
% 
% 
% B=sym(zeros(6,2));
% B(1)=0;
% B(2,:)=(jacobian(WIP.ddtheta,[u_l,u_r]));
% B(4,:)=(jacobian(WIP.ddpsi_l,[u_l,u_r]));
% B(6,:)=(jacobian(WIP.ddpsi_l,[u_l u_r]));
% B(8,:)=(jacobian(WIP.ddx,[u_l,u_r]));
% B(10,:)=(jacobian(WIP.ddy,[u_l,u_r]));
% B(12,:)=(jacobian(WIP.ddphi,[u_l u_r]));










B_=subs(Bg,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi u_l u_r],[state.',0,0]);
A_=subs(Ag,[theta dtheta psi_r dpsi_r psi_l dpsi_l x dx y dy phi dphi u_l u_r],[state.',0,0]);





A_=[[A_(1,1),A_(1,2),A_(1,11),A_(1,12)];
    [A_(2,1),A_(2,2),A_(2,11),A_(2,12)];
    [A_(11,1),A_(11,2),A_(11,11),A_(11,12)];
    [A_(12,1),A_(12,2),A_(12,11),A_(12,12)]];

B_=[B_(1,1:2);B_(2,1:2);B_(11,1:2);B_(12,1:2)];

% A=[[A(1,1),A(1,2),A(1,7),A(1,8),A(1,9),A(1,10)];
%     [A(2,1),A(2,2),A(2,7),A(2,8),A(2,9),A(2,10)];
%     [A(7,1),A(7,2),A(7,7),A(7,8),A(7,9),A(7,10)];
%     [A(8,1),A(8,2),A(8,7),A(8,8),A(8,9),A(8,10)];
%     [A(9,1),A(9,2),A(9,7),A(9,8),A(9,9),A(9,10)];
%     [A(10,1),A(10,2),A(10,7),A(10,8),A(10,9),A(10,10)]];
% 
% B=[B(1,1:2);B(2,1:2);B(7,1:2);B(8,1:2);B(9,1:2);B(10,1:2)];


end

function K = lqr__fun(state)
%A_step,B_step
[A_step,B_step]=linWip(state);
d=size(B_step)*[1;0];
horizon = 10;
P_f=eye(d);
%state=state(:).';
%state_d=[0,0,0,0,0,0];
u_l = zeros(1,horizon);
u_r = zeros(1,horizon);

%state_vec = repmat(state',1,horizon);

Q = eye(d)*1;
Q(1,1) = 5;
R = eye(2)*0.1;
P_vec = zeros(d,d*(horizon+1));

P_vec(:,d*(horizon+1) - (d-1):d*(horizon+1)) = P_f;


P_vec(:,d*(horizon+1) - (d-1):d*(horizon+1)) = P_f;
A_vec = [];
B_vec = [];
    for step = 1:horizon-1
       %linearizzo
       %state_actual = state_vec(:,horizon-step + 1);
      % [A_step,B_step] = linearization_discretization_fun(u_l(:,horizon-step),u_r(:,horizon-step),state_actual(1),state_actual(2),state_actual(3),state_actual(4));
      
       A_vec = [A_vec,A_step];
       B_vec = [B_vec,B_step];
       
       %calculate P_step
       P_next = P_vec(:,d*(horizon+2 - step) - (d-1):d*(horizon+2-step));
disp(P_next)
       Q_uu = R + B_step'*P_next*B_step; 
       P_step = Q + A_step'*P_next*A_step - (-pinv(Q_uu)*B_step'*P_next*A_step)'*Q_uu*(-pinv(Q_uu)*B_step'*P_next*A_step);
       P_vec(:,d*(horizon+1 - step) - (d-1):d*(horizon+1 - step)) = P_step;

    end
    K = (pinv(R + B_step'*P_next*B_step)*B_step'*P_next*A_step);
%     u = -K*(state()' - state_d()');
    
%     u_l = u(1);
%     u_r = u(2);
    
end