
syms F T theta DDx DDtheta Dtheta x Dx mp l mc g

eq(1)= F-T*sin(theta)==mc*DDx;
eq(2)= T*sin(theta)==mp*DDx-mp*l*DDtheta*cos(theta)+mp*l*Dtheta^2*sin(theta);
eq(3)= -mp*g-T*cos(theta)==-mp*l*DDtheta*sin(theta)-mp*l*Dtheta^2*cos(theta);

sol=(solve(eq,[DDx DDtheta T]));
disp(sol)
sol.DDx=simplify(sol.DDx);
sol.DDtheta=simplify(sol.DDtheta);
sol.T=simplify(sol.T);
disp(sol)
