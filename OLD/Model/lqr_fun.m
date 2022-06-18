function K = lqr_fun(A_step,B_step)
d=size(B_step)*[1;0];
horizon = 3;
P_f=eye(d);
%state=state(:).';
%state_d=[0,0,0,0,0,0];
u_l = zeros(1,horizon);
u_r = zeros(1,horizon);

%state_vec = repmat(state',1,horizon);

Q = eye(d)*1;
Q(1,1) = 5;
R = eye(2)*0.1;
P_vec = zeros(d,d*(horizon+1));

P_vec(:,d*(horizon+1) - (d-1):d*(horizon+1)) = P_f;


P_vec(:,d*(horizon+1) - (d-1):d*(horizon+1)) = P_f;
A_vec = [];
B_vec = [];
    for step = 1:horizon-1
       %linearizzo
       %state_actual = state_vec(:,horizon-step + 1);
      % [A_step,B_step] = linearization_discretization_fun(u_l(:,horizon-step),u_r(:,horizon-step),state_actual(1),state_actual(2),state_actual(3),state_actual(4));
      
       A_vec = [A_vec,A_step];
       B_vec = [B_vec,B_step];
       
       %calculate P_step
       P_next = P_vec(:,d*(horizon+2 - step) - (d-1):d*(horizon+2-step));
disp(P_next)
       Q_uu = R + B_step'*P_next*B_step; 
       P_step = Q + A_step'*P_next*A_step - (-pinv(Q_uu)*B_step'*P_next*A_step)'*Q_uu*(-pinv(Q_uu)*B_step'*P_next*A_step);
       P_vec(:,d*(horizon+1 - step) - (d-1):d*(horizon+1 - step)) = P_step;

    end
    K = (pinv(R + B_step'*P_next*B_step)*B_step'*P_next*A_step);
%     u = -K*(state()' - state_d()');
    
%     u_l = u(1);
%     u_r = u(2);
    
end