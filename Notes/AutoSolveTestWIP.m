syms tau Tr r mw DDpsi DDtheta theta psi g mp l Te Dtheta

%eq(1)= tau/r+Te*sin(theta)==mw*DDpsi/r;
eq(1)= tau+Tr==r^2*mw*DDpsi; 
eq(2)= -Tr+sin(theta)*l*mp*g==mp*l^2/12*DDtheta;
eq(3)= -Te*sin(theta)==mp*(DDpsi/r+l*DDtheta*cos(theta)-l*Dtheta^2*sin(theta));
eq(4)= -Te*cos(theta)-mp*g==mp*(-l*DDtheta*sin(theta)-l*Dtheta^2*cos(theta));


sol=solve(eq,[DDtheta DDpsi Te Tr]);
sol.DDtheta=simplify(sol.DDtheta);
sol.DDpsi=simplify(sol.DDpsi);
sol.Te=simplify(sol.Te);
sol.Tr=simplify(sol.Tr);
%sol.tau=simplify(sol.tau);
disp(sol)





