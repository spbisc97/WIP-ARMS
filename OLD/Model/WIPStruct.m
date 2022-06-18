function WIP=WIPStruct(Sym);
if Sym 
    syms mp 
    I_w
    WIP.mp=mp;
    WIP.w_r=I_w
end
WIP.DDtheta=(u_l + u_r + g*l*m_b*sin(theta))/(m_b*l^2 + I_b);
WIP.DDphi_r=u_r/I_w;
WIP.DDphi_l=u_l/I_w;
WIP.DDx= dpsi_r*(dpsi_r*((w_r^2*sin((w_r*(psi_l - psi_r))/w_dist))/w_dist - (w_r^3*cos((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) + (dtheta*w_r^2*sin((w_r*(psi_l - psi_r))/w_dist))/w_dist + (dpsi_l*w_r^3*cos((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) - dpsi_l*(dpsi_l*((w_r^2*sin((w_r*(psi_l - psi_r))/w_dist))/w_dist + (w_r^3*cos((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) + (dtheta*w_r^2*sin((w_r*(psi_l - psi_r))/w_dist))/w_dist - (dpsi_r*w_r^3*cos((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) + ddpsi_l*(0.5000*w_r*cos((w_r*(psi_l - psi_r))/w_dist) - (w_r^2*sin((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist) + ddpsi_r*(0.5000*w_r*cos((w_r*(psi_l - psi_r))/w_dist) + (w_r^2*sin((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist) + ddtheta*w_r*cos((w_r*(psi_l - psi_r))/w_dist) - (dtheta*w_r^2*sin((w_r*(psi_l - psi_r))/w_dist)*(dpsi_l - dpsi_r))/w_dist;
WIP.DDy= dpsi_r*(dpsi_r*((w_r^2*cos((w_r*(psi_l - psi_r))/w_dist))/w_dist + (w_r^3*sin((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) + (dtheta*w_r^2*cos((w_r*(psi_l - psi_r))/w_dist))/w_dist - (dpsi_l*w_r^3*sin((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) - dpsi_l*(dpsi_l*((w_r^2*cos((w_r*(psi_l - psi_r))/w_dist))/w_dist - (w_r^3*sin((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) + (dtheta*w_r^2*cos((w_r*(psi_l - psi_r))/w_dist))/w_dist + (dpsi_r*w_r^3*sin((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist^2) - ddpsi_l*(0.5000*w_r*sin((w_r*(psi_l - psi_r))/w_dist) + (w_r^2*cos((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist) - ddpsi_r*(0.5000*w_r*sin((w_r*(psi_l - psi_r))/w_dist) - (w_r^2*cos((w_r*(psi_l - psi_r))/w_dist)*(0.5000*psi_l + 0.5000*psi_r + theta))/w_dist) - ddtheta*w_r*sin((w_r*(psi_l - psi_r))/w_dist) - (dtheta*w_r^2*cos((w_r*(psi_l - psi_r))/w_dist)*(dpsi_l - dpsi_r))/w_dist;
WIP.DDphi= -(w_r*(ddpsi_l - ddpsi_r))/w_dist;



end