function CartPend_main_compare(Q, R, Qn, wclose)
    %if arg are less then 3 set wclose(close windows) to false
    if nargin < 4 ||isempty(wclose)
        wclose = 1;
    end
    if nargin < 2 || isempty(R)
        R=0.001;
    end
    
    if nargin <1 || isempty(Q)
        Q=diag([10,1,10,1])*0.01;
        
    end
    if nargin <3|| isempty(Q)
        Qn=diag([10,1,10,1])*100;
    end

    if (wclose)
        close all
    end

    clc;
    % global u;
    u = 0;
    y = [0; 0; pi; 0]; %initial point
    t = 0.01; %initial time
    tf = 15; %final time
    dt = 0.01; %increasing time %time step

    time = t:dt:tf+dt; %time array
    %traj_d=-0.06*(time-1.7).^5+0.6; %desired trajectory
    traj_d = repmat([1; 0; pi; 0], [1, (tf / dt)+1]);
    %traj_d(1,:)=-2*sin(time/2); %desired trajectory
    state_array = [y]; %array degli stati
    control_array = []; %array del controllo
    time_array = [t]; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = traj_d(:,1); %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")

    il=iLQR_GNMS(CartPend(),Q,R,Qn);
    il.order=[1,3,nan,nan;2,4,nan,nan];
    if coder.target("MATLAB")
    il.names=["x", "dx", "phi", "dphi"];   
    il.plot_steps=1000000;  
    il.plot_start=false;
    il.plot_end=true;
    else  
    il.names=[];   
    il.plot_steps=inf;  
    il.plot_start=false;
    il.plot_end=false;
    end
    il.plot_duration=0;
    il.defects_max=1e-5;
    il.horizon=3;

    il_ss=il;
    il_ms=il;
    il_ms_1=il;
    il_ms_1.pieces=1;


    %define plot location
    if coder.target("MATLAB")
    il_ss.plot_figure=figure("name","SS",'units','normalized','OuterPosition',[0 0  .33 1]);
    il_ms.plot_figure=figure("name","MS",'units','normalized','OuterPosition',[0.33 0  .33 1]);
    il_ms_1.plot_figure=figure("name","MS_1",'units','normalized','OuterPosition',[0.66 0  .33 1]);
    end
    il_ms.pieces=8;
    tic

    u_ms = il_ms.MS_iLQR(y,traj_d(:,1:end),t,u);
    ms_time=toc;
    tic
    u_ss = il_ss.SS_iLQR(y,traj_d(:,1:end),t,u);
    ss_time=toc;
    tic
    u_ms_1 = il_ms_1.MS_iLQR(y,traj_d(:,1:end),t,u);
    ms_1_time=toc;
    disp("times")

    disp(ms_time)
    disp(ss_time)
    disp(ms_1_time)


end
