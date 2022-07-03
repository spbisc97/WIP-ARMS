function ret=main()
    %#codegen
    ret=0;

    R=0.0001;
    Q=diag([10,1,10,1]);
    y = [-1; 0; pi; 0]; %initial point
    tf=10;
    t = 0.01; %initial time
    dt = 0.01; %increasing time %time step

    traj_d = repmat([1; 0; pi; 0], 1, (tf / dt)+1);
    traj_d(1,:) =sin(1:(tf / dt)+1);


    il=iLQR_GNMS(CartPend(),Q,R,Q);
    il.order=[1,3,nan,nan;2,4,nan,nan];
    il.defects_max=1e-3;
    il.horizon=10;

    tic
    il.MS_iLQR(y,traj_d(:,1:end),t,0);
    toc
    


    il_ms=il;
    il_ms.pieces=16;
    state_array=ones(4,il_ms.horizon_disc).*nan;
    time_array=ones(1,il_ms.horizon_disc).*nan;
    state_array(:,1) = y;
    time_array(1)=t;
    u=zeros(1,il_ms.horizon_disc-1);
    for n=1:il_ms.horizon_disc-300
        if mod(n-1,60)==0
        tic
            u(:,n:end)=il_ms.MS_iLQR(y,traj_d(:,n:end),t,u(:,n:end));
        toc
        end
        y = il_ms.dynamics_rk4(y, u(:,n));
        state_array(:,n+1) = y;
        t = t + dt; %time increment
        time_array(n+1)=t;
        plot(time_array,state_array)
        drawnow

       
    end
return