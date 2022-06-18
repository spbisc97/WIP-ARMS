clear all; 
clc; close all ;

Iw=mw*r*r;


v_w=r*Dtheta;
v_p^2=


L=1/2*mw*v_w^2+1/2*Iw*omega_w^2+1/2*mp*v_p^2 - mp*g*l*cos(theta);%L=T-U;

L=


Q=[0,0,0];
LG=EulerLagrange(vars,dvars,L,Q,0);

LG=simplify(LG);
disp(LG)

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
