%WIP with ARMS
%configuration Space
%all vector are vertical by convention
clear all%

symbolic=false;

syms x y phi %se(2) position
syms theta psi_r psi_l %body position and wheel position

syms dx dy dphi dtheta dpsi_r dpsi_l 
syms ddx ddy ddphi ddtheta ddpsi_r ddpsi_l 

syms g u_r u_l



vars= [x y phi theta psi_r psi_l ];
dvars= [dx dy dphi dtheta dpsi_r dpsi_l ];
if symbolic
%parameters 
%took similar to WIP paper with small wip (19.5cm * 10.3cm)
syms l%=5e-2; %[m] distance from wheel axis to body center of mass 
syms m_b%=3e-1; %[kg]
syms m_w%=3e-2;%[kg]
syms w_r   %=3e-2; %wheel radius
w_dm=2*w_r; %=2*3e-2; %wheel diameter

%=2*3e-2; %wheel diameter
syms w_dist%=10e-3;  %[m]distance btw wheels
%inertia
syms  I_bxx I_byy I_bzz 
syms I_wxx I_wyy I_wzz 
syms I_axx I_ayy I_azz 

else 
l=5e-2; %[m] distance from wheel axis to body center of mass 
m_b=30e-2; %[kg]
m_w=3e-2;%[kg]
w_r=3e-2; %wheel radius
w_dm=w_r; %=2*3e-2; %wheel diameter

w_dist=10e-2;  %[m]distance btw wheels
%inertia
I_bxx=550e-6 ;I_byy=8e-6; I_bzz=5e-6; 
I_wxx=5e-6; I_wyy=8e-6; I_wzz=5e-6;
syms I_axx I_ayy I_azz 
 
end 


%inertia
I_b=diag([I_bxx ,I_byy ,I_bzz]);
I_a=diag([I_axx, I_ayy, I_azz]);
I_w=diag([I_wxx, I_wyy, I_wzz]);

%body position
Bp=[x+l*sin(theta)*cos(phi);...
    y+l*sin(theta)*sin(phi);...
    (w_dm/2)+l*cos(theta)];
V_b=diff_fun(Bp,vars,dvars);%jacobian(Bp,vars)*dvars
omega_b=[-dphi*sin(theta);dtheta;dphi*cos(theta)];
%dphi has a minus bc when theta is 90 the robot lies on the xref-yref plan
%and the xbody has opposite direction to zref 


%wheel position
Wp_r=[x+w_dist/2*sin(phi);y-w_dist/2*cos(phi);w_dm/2];
Wp_l=[x-w_dist/2*sin(phi);y+w_dist/2*cos(phi);w_dm/2];
V_wl=diff_fun(Wp_l,vars,dvars);
V_wr=diff_fun(Wp_r,vars,dvars);
omega_wl=[0,dpsi_l,dphi].';% on Ixx the rotation is always zero bc of the wheel always stands up
omega_wr=[0,dpsi_r,dphi].';


%% Lagrangian
%Body K
Tb=(1/2)*(m_b*(V_b.')*V_b+omega_b.'*I_b*omega_b);




%Wheels k
Tw=(1/2)*(m_w*(V_wl.')*V_wl+omega_wl.'*I_w*omega_wl)...
    +(1/2)*(m_w*(V_wr.')*V_wr+omega_wr.'*I_w*omega_wr);


%Body P
Vb=m_b*g*Bp(3);



%Lagrangian
T=Tb+Tw;
V=Vb;

L=T-V;
disp(L)

%Sostituiamo variabili Tangenti (vincolate) se(2)
sub_phi=w_r/w_dist*(psi_r-psi_l);
sub_dphi=w_r/w_dist*(dpsi_r-dpsi_l);
% sub_x=simplify((psi_l+psi_r+2*theta)/2*w_r*cos(sub_phi));
% sub_y=simplify((psi_l+psi_r+2*theta)/2*w_r*sin(sub_phi));
sub_dx=simplify((dpsi_l+dpsi_r+2*dtheta)/2*w_r*cos(sub_phi));
sub_dy=simplify((dpsi_l+dpsi_r+2*dtheta)/2*w_r*sin(sub_phi));
L=simplify(subs(L,[dx dy dphi phi],[sub_dx sub_dy sub_dphi sub_phi]));

vars= [theta psi_r psi_l ].';
dvars= [dtheta dpsi_r dpsi_l ].';
ddvars=[ddtheta ddpsi_r ddpsi_l ].';


Q=[u_r+u_l,u_r,u_l];
LG=EulerLagrange(vars,dvars,L,Q,2);







LG=simplify(LG);
disp(LG)



disp("Other")
sub_ddx=simplify(jacobian(sub_dx,[psi_l,psi_r,theta,dpsi_l,dpsi_r,dtheta])*[dpsi_l,dpsi_r,dtheta,ddpsi_l,ddpsi_r,ddtheta].');
sub_ddy=simplify(jacobian(sub_dy,[psi_l,psi_r,theta,dpsi_l,dpsi_r,dtheta])*[dpsi_l,dpsi_r,dtheta,ddpsi_l,ddpsi_r,ddtheta].');
sub_ddphi=simplify(jacobian(sub_dphi,[psi_l,psi_r,theta,dpsi_l,dpsi_r,dtheta])*[dpsi_l,dpsi_r,dtheta,ddpsi_l,ddpsi_r,ddtheta].');
disp(sub_ddx)
disp(sub_ddx)
disp(sub_ddphi)

EQs=ddvars(:)==LG;
WIP=solve(EQs,ddvars);
WIP.ddtheta=simplify(WIP.ddtheta);
WIP.ddpsi_l=simplify(WIP.ddpsi_l);
WIP.ddpsi_r=simplify(WIP.ddpsi_r);


%put into struct also the other vars

WIP.ddx=subs(sub_ddx,[ddtheta,ddpsi_l,ddpsi_r],[WIP.ddtheta,WIP.ddpsi_l,WIP.ddpsi_r]);
WIP.ddy=subs(sub_ddy,[ddtheta,ddpsi_l,ddpsi_r],[WIP.ddtheta,WIP.ddpsi_l,WIP.ddpsi_r]);
WIP.ddphi=subs(sub_ddphi,[ddtheta,ddpsi_l,ddpsi_r],[WIP.ddtheta,WIP.ddpsi_l,WIP.ddpsi_r]);

%writestruct(WIP,"WIP.xml")
disp(WIP)


function diffun=diff_fun(fun,vars,dvars)
diffun=zeros(length(fun),1);
for i=1:length(vars)
    diffun=diffun+diff(fun,vars(i))*dvars(i);
end
end





function EQ = EulerLagrange(s,ds,L,Q,varargin)
% EQ = EulerLagrange(s,ds,L,Q,varargin) computes the Euler-Lagrange expression
% for a system, whose dynamics are defined by the Euler-Lagrange equation:
%
%     d   dL       dL
% Q = -- ----  -  ----
%     dt d(ds)     ds
%
% s   the state vector (symbols - not equations)
% ds  the state derivative vector (symbols - not equations)
% L   the lagrangian (L = E_trans + E_rotational - E_potential)
% Q   the external forces
% EQ  is the solution to the Euler-Lagrange equation. EQ is a nx1 vector
%     corresponding to the number of states in the system.
%
% The last input is a verbosity variable to silence the fprintf outputs inside
% the function. Possible inputs:
% 0   print nothing (DEFAULT),
% 1   print result,
% 2   print result and derivative terms
%

verbosity = int8(0);
if ~isempty(varargin)
    % If an extra parameter is specified, set verbosity variable
    verbosity = int8(varargin{1});
    
    % Initialize string array for holding syms expressions
    if verbosity > 0
        print_expr = strings(length(s),3);
    end
end

% Return symbolic zeros if dimensions doesnt match
if length(s) ~= length(ds)
    fprintf("<strong>[Euler-Lagrange Error]</strong> s and ds must have equal lengths! (length(s) = %g and length(ds) = %g)\n", [length(s),length(ds)]);
    return;
elseif length(s) ~= length(Q)
    fprintf("<strong>[Euler-Lagrange Error]</strong> s and Q must have equal lengths! (length(s) = %g and length(Q) = %g)\n", [length(s),length(Q)]);
    return;
end

% Ensure generalized coordinates as coloumn vectors
s = s(:);
ds = ds(:);

% Make variables for higher derivatives
dds = str2sym( "d" + string(ds) );
if verbosity > 1
    fprintf("Noting higher derivatives as: ")
    fprintf("dds = [ "); fprintf('%s ', string(dds)); fprintf("];\n");
end

% Compute general equations, EQ, for all generalized coordinates, s:
EQ = sym(zeros(length(s),1));
for ii = 1:length(s)
    
    % Partial Derivatives
    partial_s  = diff(L, s(ii));
    partial_ds = diff(L,ds(ii));

    % Time derivatives (applying the chain-rule df(x,t)/dt = df/dx * dx/dt)
    partial_dt_ds = jacobian(partial_ds, [s;ds]) * [ds;dds];
    
    % Solve Euler-Lagrange equation with input forces, Q
    EQ(ii) = reduce( solve(Q(ii) == partial_dt_ds - partial_s, dds(ii)) );
    if isempty(EQ(ii))
        fprintf("<strong>[Euler-Lagrange Warning]</strong> State %i did not have a nonzero solution.\n", ii);
    end
    
    % Save equations for printing
    if verbosity > 0
        if verbosity > 1
            print_expr(ii,1) = sprintf( ['%g. ',char(reduce(partial_s))], ii);
            print_expr(ii,2) = sprintf( ['%g. ',char(reduce(partial_dt_ds))], ii);
        end
        print_expr(ii,3) = sprintf( ['%g. ',char(dds(ii)),' = ',char(EQ(ii))], ii);
    end
end

% Print expressions
if verbosity > 0
    if verbosity > 1
        fprintf("------------------------------------\nDerivative(s) of the potential term:\n------------------------------------\n")
        fprintf(1, '%s \n', print_expr(:,1)); fprintf("\n");
        
        fprintf("------------------------------------\nDerivative(s) of the kinetic term:\n------------------------------------\n")
        fprintf(1, '%s \n', print_expr(:,2)); fprintf("\n");
    end
    fprintf("------------------------------------\nGeneral Equation(s):\n------------------------------------\n")
    fprintf(1, '%s \n', print_expr(:,3)); fprintf("\n");
end

end

% Utility function
function expr = reduce(expr)
expr = simplify(expand(expr));
end
