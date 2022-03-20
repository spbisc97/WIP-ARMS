syms mw DDpsi DDtheta psi theta Dtheta Dpsi l mp r g tau
eq(1)= mw*DDpsi*r==tau/r + mp*(cos(theta)*DDtheta+sin(theta)*Dtheta^2)*l;
eq(2)=(mp*l^2+mp*l^2/12)*DDtheta==mp*g*sin(theta)*l-mp*DDtheta*r*cos(theta)*l;
sol=solve(eq,[DDpsi ;DDtheta]);
