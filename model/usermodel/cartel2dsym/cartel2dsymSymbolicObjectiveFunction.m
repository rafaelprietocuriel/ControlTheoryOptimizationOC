function out=cartel2dsymSymbolicObjectiveFunction(arcid,s)
% returns the objective function
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
out=mystr2sym('exp(-r*t)*(-w*s*C1*C2-(1-w)*h*(C1+C2)-gamma1*u1-gamma2*u2-kappa*v)');
	
if s
	ctrl=cartel2dsymSymbolicOptimalControl(arcid);
	[lagmcc,lagmsc]=cartel2dsymSymbolicLagrangeMultiplier(arcid);
	out=subs(out,{'u1','u2','v','lagmcc1','lagmcc2','lagmcc3','lagmcc4'},{ctrl(1),ctrl(2),ctrl(3),lagmcc(1),lagmcc(2),lagmcc(3),lagmcc(4)});
end
