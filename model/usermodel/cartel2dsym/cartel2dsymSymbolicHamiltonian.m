function out=cartel2dsymSymbolicHamiltonian(arcid,s)
% returns the Hamiltonian value for cartel2dsym
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
out=mystr2sym('-w*s*C1*C2-(1-w)*h*(C1+C2)-gamma1*u1-gamma2*u2-kappa*v+lambda1*(tau+(b+(delta-b)*exp(-xi*v))*rho*C1-eta*u1^alpha*(C1+epsilon)^beta-theta*C1*C2-omega*C1^2)+lambda2*(tau+(b+(delta-b)*exp(-xi*v))*rho*C2-eta*u2^alpha*(C2+epsilon)^beta-theta*C1*C2-omega*C2^2)+lagmcc1*(u1-(ulow))+lagmcc2*(u2-(ulow))+lagmcc3*(v-(vlow))+lagmcc4*(B-(u1+u2+v))');
	
if s
	ctrl=cartel2dsymSymbolicOptimalControl(arcid);
	[lagmcc,lagmsc]=cartel2dsymSymbolicLagrangeMultiplier(arcid);
	out=subs(out,{'u1','u2','v','lagmcc1','lagmcc2','lagmcc3','lagmcc4'},{ctrl(1),ctrl(2),ctrl(3),lagmcc(1),lagmcc(2),lagmcc(3),lagmcc(4)});
end
