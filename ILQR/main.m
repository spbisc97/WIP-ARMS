function n=main()
    %#codegen

    R=0.0001;
    Q=diag([10,1,10,1]);
    u = 0;
    y = [-1; 0; pi; 0]; %initial point
    t = 0.01; %initial time
    tf = 15; %final time
    dt = 0.01; %increasing time %time step

    traj_d = repmat([1; 0; pi; 0], [1, (tf / dt)+1]);
    traj_d(1,:) =sin(1:(tf / dt)+1);


    il=iLQR_GNMS(CartPend(),Q,R,Q);
    il.order=[1,3,nan,nan;2,4,nan,nan];
    il.defects_max=1e-5;
    il.horizon=10;

    il_ss=il;
    il_ms=il;
    il_ms_1=il;
    il_ms_1.pieces=1;



    %il.pieces=16;
    tic;

    il_ms.MS_iLQR(y,traj_d(:,1:end),t,u);
    toc;
    tic;
    il_ss.SS_iLQR(y,traj_d(:,1:end),t,u);
    toc;
    tic;
    il_ms_1.MS_iLQR(y,traj_d(:,1:end),t,u);
    toc;

    n=1;
   