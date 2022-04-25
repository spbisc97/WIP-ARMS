function u=LQR_function(state,state_d,Q,R)

u=0;

 
horizon=100;
P_f=eye(1);
P=P_f;
%state_d=state_d-state;

for step = 1:horizon-1
    [A,B]=linearization_discretization(u,state);
    P_next=A'*P*A-(A'*P*B)*pinv(R+B'*P*B)*(B'*P*A)+Q; 
    P=P_next;
    
end
K=pinv(R+B'*P*B)*(B'*P*A);
disp('')


u=-K*(state-state_d);

end
