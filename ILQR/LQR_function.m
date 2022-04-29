function u=LQR_function(state,state_d,Q,R)

u=0;

 
horizon=10;
P_f=eye(1);
P=P_f;

for step = 1:horizon-1
    [A,B]=linearization_discretization(u,state);%linearizziamo il sistema nel punto 
    P_next=A'*P*A-(A'*P*B)*pinv(R+B'*P*B)*(B'*P*A)+Q;  %discrete ARE solution 
    P=P_next;%next P for riccati 
end
K=pinv(R+B'*P*B)*(B'*P*A);
disp('')


u=K*(state_d-state);

end
