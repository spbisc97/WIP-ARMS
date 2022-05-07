function [A, B] = linearization_discretization(u, state)
    y = state;
    %y=[0 0 pi 0]; choose linearization point different from state
    y1 = y(1); %posizione x
    y2 = y(2); %velocità x
    y3 = y(3); %angolo phi
    y4 = y(4); %velocità angolare phi
    Ts = 0.01;
    m = 1;
    M = 5;
    L = 2;
    g = -10;
    d = 1;

    % sys=[y(2),(1/D)*(-m^2*L^2*g*cos(y(3))*sin(y(3)) + m*L^2*(m*L*y(4)^2*sin(y(3)) - d*y(2))) + m*L*L*(1/D)*u,y(4),(1/D)*((m+M)*m*g*L*sin(y(3)) - m*L*cos(y(3))*(m*L*y(4)^2*sin(y(3)) - d*y(2))) - m*L*cos(y(3))*(1/D)*u]

    A = [
        [0, 1, 0, 0];
        [0, -d / (M - m * (cos(y3)^2 - 1)), (m * sin(2 * y3) * (- 2 * L * m * sin(y3) * y4^2 + 2 * d * y2 + g * m * sin(2 * y3))) / (2 * (M + m / 2 - (m * cos(2 * y3)) / 2)^2) - (4 * m * u * sin(2 * y3)) / (2 * M + m - m * cos(2 * y3))^2 - (m * (- L * cos(y3) * y4^2 + g * cos(2 * y3))) / (M + m / 2 - (m * cos(2 * y3)) / 2), (2 * L * m * y4 * sin(y3)) / (m * sin(y3)^2 + M)];
        [0, 0, 0, 1];
        [0, (d * cos(y3)) / (L * (M - m * (cos(y3)^2 - 1))), (u * sin(y3)) / (L * (M - m * (cos(y3)^2 - 1))) - (L * m * sin(y3) * (- L * m * sin(y3) * y4^2 + d * y2) + L^2 * m^2 * y4^2 * cos(y3)^2 - L * g * m * cos(y3) * (M + m)) / (L^2 * m * (M - m * (cos(y3)^2 - 1))) - (2 * cos(y3) * sin(y3) * (L * m * cos(y3) * (- L * m * sin(y3) * y4^2 + d * y2) + L * g * m * sin(y3) * (M + m))) / (L^2 * (M - m * (cos(y3)^2 - 1))^2) + (2 * m * u * cos(y3)^2 * sin(y3)) / (L * (M - m * (cos(y3)^2 - 1))^2), - (2 * m * y4 * cos(y3) * sin(y3)) / (M - m * (cos(y3)^2 - 1))]; ];

    B = [0
        1 / (M - m * (cos(y3)^2 - 1))
        0
        - cos(y3) / (L * (M - m * (cos(y3)^2 - 1)))];
    n=4;
    B = pinv(eye(n) - A * Ts / 2.) * B * sqrt(Ts);
    A = (eye(n) + A * Ts / 2.) * pinv(eye(n) - A * Ts / 2.);%tustin, bilinear trans

end
