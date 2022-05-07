function new_state = euler_integration_fun(state, dy, dt)
    %integrate as euler
    %%test with single value
    %new_state=state+dy*dt;
    %%test for inverted pendulum
    new_state(2) = state(2) + dy(2) * dt;
    new_state(1) = state(1) + new_state(2) * dt;

    new_state(4) = state(4) + dy(4) * dt;
    new_state(3) = state(3) + new_state(4) * dt;
    new_state=new_state.';
end
