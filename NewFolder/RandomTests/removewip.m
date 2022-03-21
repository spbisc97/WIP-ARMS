function dy = removewip(y,mp,mc,l,g,d,r,F)
theta=y(3);
Dtheta=y(4);
mw=mc;
tau=F;


dy(1,1) = y(2);
dy(2,1) = (r*(tau - g*mp*sin(theta) + g*l*mp*sin(theta)))/(mw*r^2 - mp*cos(theta));
dy(3,1) = y(4);
dy(4,1) = -(- g*mw*sin(theta)*r^2 + tau*cos(theta) + g*l*mp*cos(theta)*sin(theta))/(l*(mw*r^2 - mp*cos(theta)));

disp(dy.')
%pause