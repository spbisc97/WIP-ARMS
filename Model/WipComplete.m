
clear all, close all, clc


syms x theta Dtheta Dx
syms DDx DDtheta  F mp mc l g 

syms Tr r mw g mp l Te  tau
syms DDpsi DDtheta theta Dtheta psi

g=9.81;  l=2;  mc=5;
mp=2; mw=1; r=1;
tau=F;
Iw=mw*r*r;
Ip=mp*l*l/12;
stab=0;
eq(1)= (mp+mw+Iw/r^2)*DDpsi*r+mp*l*cos(theta)*DDtheta-mp*l*sin(theta)*Dtheta^2==tau/r; 
eq(2)= (mp*l*cos(theta))*DDpsi*r+(Ip+mp*l)*DDtheta-mp*l*g*sin(theta)==-tau;







sol=solve(eq,[DDtheta DDpsi ]); %%add as much var as needed
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

Q = [10 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 100];
R = 70;

%%
det(ctrb(A,B))

%%
K =lqr(A,B,Q,R);
s=-1;
disp("start Int")
tspan = 0:.001:100;
if(s==-1)
    y0 = [0; 0; stab+0.1; 0];
    [t,y] = ode45(@(t,y)removewip(y,mp,mc,l,g,0,r,-K*(y-[0; 0; stab; 0])),tspan,y0);
elseif(s==1)
    y0 = [-3; 0; 0; 0];
% % [t,y] = ode45(@(t,y)((A-B*K)*(y-[0; 0; pi; 0])),tspan,y0);
    [t,y] = ode45(@(t,y)removewip(y,mp,mc,l,g,0,r,-K*(y-[0; 0; 0; 0])),tspan,y0);
else
    
end

for k=1:100:length(t)
    figure(1)
    drawwip(y(k,:),mp,mc,l,r);
end

%% 

function drawwip(y,m,~,L,r)
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

function dy = removewip(y,mp,mw,l,g,d,r,F)
theta=y(3);
Dtheta=y(4);
tau=F;


dy(1,1) = y(2);
dy(2,1) = (12*tau + l*tau + 12*r*tau*cos(theta) + 12*Dtheta^2*l*mp*r*sin(theta) + Dtheta^2*l^2*mp*r*sin(theta) - 12*g*l*mp*r*cos(theta)*sin(theta))/(r^2*(12*mp + 24*mw + l*mp + 2*l*mw - 12*l*mp*cos(theta)^2));
dy(3,1) = y(4);
dy(4,1) = -(12*(mp*r*tau + 2*mw*r*tau + l*mp*tau*cos(theta) - g*l*mp^2*r*sin(theta) - 2*g*l*mp*mw*r*sin(theta) + Dtheta^2*l^2*mp^2*r*cos(theta)*sin(theta)))/(l*mp*r*(12*mp + 24*mw + l*mp + 2*l*mw - 12*l*mp*cos(theta)^2));

disp(dy.')
% pause
end