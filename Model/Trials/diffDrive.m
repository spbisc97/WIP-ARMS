clear all
clc
d=0.5 %[m]

syms x dx y dy phi dphi
state=[x,dx,y,dy,phi,dphi];%,psi_l,dpsi_l,psi_r,dpsi_r





function dy=(state,u)
    T_l=u(2);
    T_r=u(1);

    ddpsi_l=I*T_l
    ddpsi_r=I*T_r
    ddphi=(ddpsi_r-ddpsi_r)/d
    ddy=
    ddx=
    

end