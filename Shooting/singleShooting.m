
function output = shooting_method(input)
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
        
end


function F = solver(x)
    options=odeset("RelTol",1e-8,'AbsTol',[1e-8,1e-8])
    [t,u]=ode45(@equation,[0 1],[1 x],options)

    s=lenght(t)
    F=u(1)
end


function dy = equation(y,t)
    dy=zeros(2,1) 
    dy(1)=y(2)
    dy(2)=t*exp(2*t)+y(1)-t*(y(2)^2)
end