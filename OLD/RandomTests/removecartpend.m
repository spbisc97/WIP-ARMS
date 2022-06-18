function dy = removecartpend(y,mp,mc,l,g,d,F)
theta=y(3);
Dtheta=y(4);


dy(1,1) = y(2);
dy(2,1) = (- l*mp*sin(theta)*Dtheta^2 + F + g*mp*cos(theta)*sin(theta))/(- mp*cos(theta)^2 + mc + mp);
dy(3,1) = y(4);
dy(4,1) = (- l*mp*cos(theta)*sin(theta)*Dtheta^2 + F*cos(theta) + g*mc*sin(theta) + g*mp*sin(theta))/(l*(- mp*cos(theta)^2 + mc + mp));
