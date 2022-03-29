function WIP= WIP3dModel(g)


syms ddpsi ddtheta dtheta theta dpsi_l dpsi_r psi_r psi_l
syms u_l u_r 


WIP.ddtheta = (5.0000e-04*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40);
WIP.ddpsi_r = -(0.0080*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);
WIP.ddpsi_l = -(0.0080*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);
WIP.ddx=(1.5000e-05*cos(0.3000*psi_l - 0.3000*psi_r)*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40) - (1.2000e-04*cos(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - (1.2000e-04*cos(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - 0.3000*dpsi_l*sin(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta) + 0.3000*dpsi_r*sin(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta);
WIP.ddy=(1.2000e-04*sin(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) + (1.2000e-04*sin(0.3000*psi_l - 0.3000*psi_r)*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - (1.5000e-05*sin(0.3000*psi_l - 0.3000*psi_r)*(1.0889e+46*u_l + 1.0889e+46*u_r + 4.9001e+42*dtheta^2*sin(theta) + 1.3484e+43*dpsi_l^2*sin(2*theta) + 1.3484e+43*dpsi_r^2*sin(2*theta) - 6.8907e+43*dtheta^2*sin(2*theta) - 3.0625e+47*u_l*cos(theta) - 3.0625e+47*u_r*cos(theta) + 3.4709e+45*g*sin(theta) - 2.6969e+43*dpsi_l*dpsi_r*sin(2*theta)))/(4.9001e+39*cos(theta) - 6.8907e+40*cos(theta)^2 + 8.9462e+40) - 0.3000*dpsi_l*cos(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta) + 0.3000*dpsi_r*cos(0.3000*psi_l - 0.3000*psi_r)*(0.0150*dpsi_l + 0.0150*dpsi_r + 0.0300*dtheta);
WIP.ddphi=(0.0024*(3.1739e+63*dpsi_l^2*sin(theta) - 5.6642e+67*u_r - 3.1032e+68*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 2.1808e+68*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 5.9193e+67*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 7.9074e+67*u_l*cos(theta) - 5.6476e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) + 1.3169e+63*dpsi_l*dtheta*sin(theta) - 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) + 2.9566e+64*dpsi_l*dtheta*sin(2*theta) + 1.3169e+63*dpsi_l*dtheta*sin(3*theta) - 9.2593e+63*dpsi_l*dtheta*sin(4*theta) - 2.9566e+64*dpsi_r*dtheta*sin(2*theta) - 1.3169e+63*dpsi_r*dtheta*sin(3*theta) + 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62) - (0.0024*(3.1739e+63*dpsi_l^2*sin(theta) - 3.1032e+68*u_r - 5.6642e+67*u_l + 3.1739e+63*dpsi_r^2*sin(theta) - 1.1376e+65*dtheta^2*sin(theta) + 5.9193e+67*u_l*cos(2*theta) + 4.1152e+67*u_l*cos(3*theta) + 2.1808e+68*u_r*cos(2*theta) + 4.1152e+67*u_r*cos(3*theta) + 8.1695e+65*g*sin(2*theta) - 4.4444e+65*g*sin(3*theta) - 3.0864e+65*g*sin(4*theta) + 4.5704e+63*dpsi_l^2*sin(2*theta) + 1.9748e+63*dpsi_l^2*sin(3*theta) - 1.7267e+63*dpsi_l^2*sin(4*theta) - 1.1991e+63*dpsi_l^2*sin(5*theta) + 4.5704e+63*dpsi_r^2*sin(2*theta) + 1.9748e+63*dpsi_r^2*sin(3*theta) - 1.7267e+63*dpsi_r^2*sin(4*theta) - 1.1991e+63*dpsi_r^2*sin(5*theta) - 2.4509e+64*dtheta^2*sin(2*theta) + 3.1193e+64*dtheta^2*sin(3*theta) + 9.2593e+63*dtheta^2*sin(4*theta) - 5.6476e+67*u_l*cos(theta) - 7.9074e+67*u_r*cos(theta) + 1.6209e+66*g*sin(theta) - 6.3477e+63*dpsi_l*dpsi_r*sin(theta) - 1.3169e+63*dpsi_l*dtheta*sin(theta) + 1.3169e+63*dpsi_r*dtheta*sin(theta) - 9.1407e+63*dpsi_l*dpsi_r*sin(2*theta) - 3.9496e+63*dpsi_l*dpsi_r*sin(3*theta) + 3.4533e+63*dpsi_l*dpsi_r*sin(4*theta) + 2.3981e+63*dpsi_l*dpsi_r*sin(5*theta) - 2.9566e+64*dpsi_l*dtheta*sin(2*theta) - 1.3169e+63*dpsi_l*dtheta*sin(3*theta) + 9.2593e+63*dpsi_l*dtheta*sin(4*theta) + 2.9566e+64*dpsi_r*dtheta*sin(2*theta) + 1.3169e+63*dpsi_r*dtheta*sin(3*theta) - 9.2593e+63*dpsi_r*dtheta*sin(4*theta)))/(4.8955e+61*cos(theta) - 1.4578e+63*cos(theta)^2 - 4.2140e+61*cos(theta)^3 + 5.9259e+62*cos(theta)^4 + 8.9379e+62);

end