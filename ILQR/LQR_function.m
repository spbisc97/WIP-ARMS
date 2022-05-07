function u = LQR_function(state, state_d, Q, R)

    u = 0;

    horizon = 200000000;%montepremi superenalotto 7/05/22
    P = Q;

    [A, B] = linearization_discretization(u, state); %linearizziamo il sistema nel punto

    for step = 1:horizon
        P_next = A' * P * A - (A' * P * B) * pinv(R + B' * P * B) * (B' * P * A) + Q; %discrete ARE solution
        if(P-P_next)<1e-4 %threshold P evolution
            %disp(step) 
            break %exit from the loop
        end
        P = P_next; %next P for riccati
    end
    K = pinv(R + B' * P * B) * (B' * P * A); %calculate K
    u = K * (state_d - state); %calc control
end
