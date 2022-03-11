%WIP with ARMS
%configuration Space
%all vector are vertical by convention
clear all%
syms phi %se(2) position
syms alpha_r alpha_l % position of arm is respect to body as robotics 1 manipulators

syms dphi  dalpha_r dalpha_l

vars= [phi alpha_r alpha_l];
dvars= [dphi  dalpha_r dalpha_l];


%parameters 
%took similar to WIP paper with small wip (19.5cm * 10.3cm)
syms l%=5e-2; %[m] distance from wheel axis to body center of mass 
syms l_a%=2e-2; %[m]lenght of arms
syms m_b%=3e-1; %[kg]
syms m_a%=3e-2;%[kg]
syms m_w%=3e-2;%[kg]
syms w_d%=2*3e-2; %wheel diameter

syms w_dist%=10e-3;  %[m]distance btw wheels
syms d_a%=10e-3;  %[m]distance btw arms
syms g

%inertia
% I_bxx=550e-6 I_byy=480e-6 I_bzz=150e-6 
syms I_b%=diag(I_bxx ,I_byy ,I_bzz);
% I_axx=m_a*l_a/3 I_ayy I_azz 
syms I_a%=diag(I_axx, I_ayy, I_azz);
% I_wxx I_wyy I_wzz 
syms I_w%=diag(I_wxx, I_wyy, I_wzz);

%body position
Bp=[0;...
    0;...
    0];
V_b=diff_fun(Bp,vars,dvars)%jacobian(Bp,vars)*dvars
omega_b=[0;0;dphi];


Ap_l=[0-d_a/2*sin(phi)+((l_a/2)*sin(alpha_l))*cos(phi);...
    0+d_a/2*cos(phi)+((l_a/2)*sin(alpha_l))*sin(phi);...
    l+(l_a/2)*cos(alpha_l)];
V_al=(diff_fun(Ap_l,vars,dvars))
omega_al=[dalpha_l,0,dphi].';           

Ap_r=[0+d_a/2*sin(phi)+((l_a/2)*sin(alpha_r))*cos(phi);...
    0-d_a/2*cos(phi)+((l_a/2)*sin(alpha_r))*sin(phi);...
    l+(l_a/2)*cos(alpha_r)];
V_ar=(diff_fun(Ap_r,vars,dvars))
omega_ar=[dalpha_r,0,dphi].';




%% Lagrangian
%Body K
Tb=(1/2)*(m_b*(V_b.')*V_b+omega_b.'*I_b*omega_b);


%Arms K
Ta=(1/2)*(m_a*(V_al.')*V_al+(omega_al.')*I_a*omega_al)...
    +(1/2)*(m_a*(V_ar.')*V_ar+(omega_ar.')*I_a*omega_ar);




%Body P
Vb=m_b*g*Bp(3);
%Arm P
Va=m_a*g*Ap_l(3)+m_a*g*Ap_r(3); 


%Lagrangian
T=Tb+Ta;
V=Va+Vb;

L=T-V;

Q=[0,0,0];
EulerLagrange(vars,dvars,L,Q,1)



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
