
function singleShooting
%     clc
%     clear all;
%     x=0.5;
%     fun=@solver;
%     fzero(fun,x)       
solver(0.5)
end


function F=solver(x)
    options=odeset('RelTol',1e-8,'AbsTol',[1e-8,1e-8]);
    [t,u]=ode45(@equation,[0 1],[1 0.5],options);
    %s=lenght(t);
    %F=u(s,1)-exp(1);
    %figure(1);
    %plot(t,u(:,1));
    F=x
end


function dy = equation(y,t)
    dy=zeros(2,1);
    dy(1)=y(2);
    dy(2)=t*exp(2*t)+y(1)-t*(y(2)^2);
end