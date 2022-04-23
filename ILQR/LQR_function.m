function u=LQR_function(state,state_d)

u=0;

Q = 100;
R = 0.01;
horizon=10;
P_f=eye(1);
P=P_f;


for step = 1:horizon-1
    disp(P)
    [A,B]=linearization_discretization(u,state);
    P_next=A'*P*A-(A'*P*B)*pinv(R+B'*P*B)*(B'*P*A)+Q; 
    P=P_next;
    
end
K=(R+B'*P*B)*(B'*P*A);
disp('')


u=-K*(state-state_d);

end
