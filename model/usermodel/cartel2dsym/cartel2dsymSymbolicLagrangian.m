function out=cartel2dsymSymbolicLagrangian(arcid)
%
% returns the symbolic extended Hamiltonian (Lagrangian) value for cartel2dsym
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
out=mystr2sym('-w*s*C1*C2-(1-w)*h*(C1+C2)-gamma1*u1-gamma2*u2-kappa*v+lambda1*(tau+(b+(delta-b)*exp(-xi*v))*rho*C1-eta*u1^alpha*(C1+epsilon)^beta-theta*C1*C2-omega*C1^2)+lambda2*(tau+(b+(delta-b)*exp(-xi*v))*rho*C2-eta*u2^alpha*(C2+epsilon)^beta-theta*C1*C2-omega*C2^2)+lagmcc1*(u1-(ulow))+lagmcc2*(u2-(ulow))+lagmcc3*(v-(vlow))+lagmcc4*(B-(u1+u2+v))');
	
[lagmcc,lagmsc]=cartel2dsymSymbolicLagrangeMultiplier(arcid);
	for ii=1:4
		if lagmcc(ii)==0
			out=subs(out,['lagmcc' num2str(ii)],0);
		end
	end
