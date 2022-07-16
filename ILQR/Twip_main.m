function Twip_main(Q, R, wclose)
    %if arg are less then 3 set wclose(close windows) to false
    if nargin < 3
        wclose = 1;
    end
    if nargin < 2
        R=diag([0.01,0.01]);
    end
    if nargin <1
        Q=diag([10,1,10,1,10,1]);
        Qn=diag([10,1,10,1,10,1])*1000;
    end

    if (wclose)
        close all
    end

    clc;
    % global u;
    u = [0,0,0;0,0,0];
    y = [0; 0; 0; 0;0;0]; %initial point
    t = 0.01; %initial time
    tf = 15; %final time
    dt = 0.01; %increasing time %time step

    time = t:dt:tf+dt; %time array
    %traj_d=-0.06*(time-1.7).^5+0.6; %desired trajectory
    traj_d = repmat([0; 0; 1; 0;0;0], [1, (tf / dt)+1]);
    traj_d(5,:)=-2*sin(time); %desired trajectory
    traj_d(6,:)=-2*cos(time); %desired trajectory

    state_array = [y]; %array degli stati
    control_array = []; %array del controllo
    time_array = [t]; %array del tempo (dovrebbe coincidere con la l'array "time")
    y_d_array = traj_d(:,1); %array della traiettoria (dovrebbe coincidere con la l'array "traj_d")
    il=iLQR_GNMS(Twip(),Q,R,Qn);
    il.order=[1,2,nan,nan,nan,nan;3,4,nan,nan,nan,nan;5,6,nan,nan,nan,nan];
    il.names=["phi","dphi","x" "dx", "theta","dtheta"];
    il.pieces=1;
    il.plot_steps=100000;
    il.plot_start=false;
    il.plot_end=true;
    il.plot_duration=0;
    il.defects_max=1e-7;

    il_ss=il;
    il_ms=il;
    il_ms_1=il;
    il_ms_1.pieces=1;


    %define plot location
    il_ss.plot_figure=figure("name","SS",'units','normalized','OuterPosition',[0 0  .33 1]);
    il_ms.plot_figure=figure("name","MS",'units','normalized','OuterPosition',[0.33 0  .33 1]);
    il_ms_1.plot_figure=figure("name","MS_1",'units','normalized','OuterPosition',[0.66 0  .33 1]);

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