function out=cartel2dsymSymbolicCanonicalSystem(arcid,s)
% returns expression for the canonical systems dynamics for cartel2dsym
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
out=mystr2sym(['[tau+(b+(delta-b)*exp(-xi*v))*rho*C1-eta*u1^alpha*(C1+epsilon)^beta-theta*C1*C2-omega*C1^2;' ...
	'tau+(b+(delta-b)*exp(-xi*v))*rho*C2-eta*u2^alpha*(C2+epsilon)^beta-theta*C1*C2-omega*C2^2;' ...
	'r*lambda1-(-w*s*C2-(1-w)*h+lambda1*((b+(delta-b)*exp(-xi*v))*rho-eta*u1^alpha*(C1+epsilon)^beta*beta/(C1+epsilon)-theta*C2-2*omega*C1)-lambda2*theta*C2);' ...
	'r*lambda2-( -w*s*C1-(1-w)*h-lambda1*theta*C1+lambda2*((b+(delta-b)*exp(-xi*v))*rho-eta*u2^alpha*(C2+epsilon)^beta*beta/(C2+epsilon)-theta*C1-2*omega*C2))]']);
	
switch arcid
	case 0
		outae=mystr2sym('-kappa-lambda1*(delta-b)*xi*exp(-xi*v)*rho*C1-lambda2*(delta-b)*xi*exp(-xi*v)*rho*C2');
	
	case 1
		outae=mystr2sym('-kappa-lambda1*(delta-b)*xi*exp(-xi*v)*rho*C1-lambda2*(delta-b)*xi*exp(-xi*v)*rho*C2');
	
	case 2
		outae=mystr2sym('-kappa-lambda1*(delta-b)*xi*exp(-xi*v)*rho*C1-lambda2*(delta-b)*xi*exp(-xi*v)*rho*C2');
	
	case 3
		outae=mystr2sym([]);
	
	case 4
		outae=mystr2sym(['[-gamma1-lambda1*eta*u1^alpha*alpha/u1*(C1+epsilon)^beta-lagmcc4;' ...
			' -gamma2-lambda2*eta*u2^alpha*alpha/u2*(C2+epsilon)^beta-lagmcc4]']);
	
	case 5
		outae=mystr2sym('-kappa-lambda1*(delta-b)*xi*exp(-xi*v)*rho*C1-lambda2*(delta-b)*xi*exp(-xi*v)*rho*C2');
	
	case 6
		outae=mystr2sym([]);
	
	case 7
		outae=mystr2sym('-gamma2-lambda2*eta*u2^alpha*alpha/u2*(C2+epsilon)^beta-lagmcc4');
	
	case 8
		outae=mystr2sym([]);
	
	case 9
		outae=mystr2sym('-gamma1-lambda1*eta*u1^alpha*alpha/u1*(C1+epsilon)^beta-lagmcc4');
	
	case 10
		outae=mystr2sym('-gamma2-lambda2*eta*u2^alpha*alpha/u2*(C2+epsilon)^beta-lagmcc4');
	
	case 11
		outae=mystr2sym([]);
	
	case 12
		outae=mystr2sym([]);
	
	case 13
		outae=mystr2sym([]);
	
	case 14
		outae=mystr2sym([]);
	
end
out=[out;outae];
	
