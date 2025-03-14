function [lagmcc,lagmsc]=cartel2dsymSymbolicLagrangeMultiplier(arcid)
% returns the symbolic expressions for the Lagrangian multipliers for inequality constraints for cartel2dsym
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
lagmcc=cartel2dsymSymbolicControlLagrangeMultiplier(arcid);
lagmsc=cartel2dsymSymbolicStateLagrangeMultiplier(arcid);
	
%-------------------------------------------------------------------------
function out=cartel2dsymSymbolicControlLagrangeMultiplier(arcid)
% returns the symbolic expressions for the Lagrangian multipliers for (mixed) control inequality constraints for cartel2dsym
	
switch arcid
	case 0
		out=mystr2sym(['[0,' ...
			'0,' ...
			'0,' ...
			'0]']);
	
	case 1
		out=mystr2sym(['[(gamma1*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta)/ulow,' ...
			'0,' ...
			'0,' ...
			'0]']);
	
	case 2
		out=mystr2sym(['[0,' ...
			'(gamma2*ulow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta)/ulow,' ...
			'0,' ...
			'0]']);
	
	case 3
		out=mystr2sym(['[0,' ...
			'0,' ...
			'kappa+lambda1*xi*exp(-xi*vlow)*rho*C1*delta-lambda1*xi*exp(-xi*vlow)*rho*C1*b+lambda2*xi*exp(-xi*vlow)*rho*C2*delta-lambda2*xi*exp(-xi*vlow)*rho*C2*b,' ...
			'0]']);
	
	case 4
		out=mystr2sym(['[0,' ...
			'0,' ...
			'0,' ...
			'-kappa-lambda1*xi*exp(-xi*(B-u1-u2))*rho*C1*delta+lambda1*xi*exp(-xi*(B-u1-u2))*rho*C1*b-lambda2*xi*exp(-xi*(B-u1-u2))*rho*C2*delta+lambda2*xi*exp(-xi*(B-u1-u2))*rho*C2*b]']);
	
	case 5
		out=mystr2sym(['[(gamma1*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta)/ulow,' ...
			'(gamma2*ulow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta)/ulow,' ...
			'0,' ...
			'0]']);
	
	case 6
		out=mystr2sym(['[(gamma1*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta)/ulow,' ...
			'0,' ...
			'kappa+lambda1*xi*exp(-xi*vlow)*rho*C1*delta-lambda1*xi*exp(-xi*vlow)*rho*C1*b+lambda2*xi*exp(-xi*vlow)*rho*C2*delta-lambda2*xi*exp(-xi*vlow)*rho*C2*b,' ...
			'0]']);
	
	case 7
		out=mystr2sym(['[(gamma1*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta-ulow*kappa-ulow*lambda1*xi*exp(-xi*(B-ulow-u2))*rho*C1*delta+ulow*lambda1*xi*exp(-xi*(B-ulow-u2))*rho*C1*b-ulow*lambda2*xi*exp(-xi*(B-ulow-u2))*rho*C2*delta+ulow*lambda2*xi*exp(-xi*(B-ulow-u2))*rho*C2*b)/ulow,' ...
			'0,' ...
			'0,' ...
			'-kappa-lambda1*xi*exp(-xi*(B-ulow-u2))*rho*C1*delta+lambda1*xi*exp(-xi*(B-ulow-u2))*rho*C1*b-lambda2*xi*exp(-xi*(B-ulow-u2))*rho*C2*delta+lambda2*xi*exp(-xi*(B-ulow-u2))*rho*C2*b]']);
	
	case 8
		out=mystr2sym(['[0,' ...
			'(gamma2*ulow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta)/ulow,' ...
			'kappa+lambda1*xi*exp(-xi*vlow)*rho*C1*delta-lambda1*xi*exp(-xi*vlow)*rho*C1*b+lambda2*xi*exp(-xi*vlow)*rho*C2*delta-lambda2*xi*exp(-xi*vlow)*rho*C2*b,' ...
			'0]']);
	
	case 9
		out=mystr2sym(['[0,' ...
			'(gamma2*ulow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta-ulow*kappa-ulow*lambda1*xi*exp(-xi*(B-u1-ulow))*rho*C1*delta+ulow*lambda1*xi*exp(-xi*(B-u1-ulow))*rho*C1*b-ulow*lambda2*xi*exp(-xi*(B-u1-ulow))*rho*C2*delta+ulow*lambda2*xi*exp(-xi*(B-u1-ulow))*rho*C2*b)/ulow,' ...
			'0,' ...
			'-kappa-lambda1*xi*exp(-xi*(B-u1-ulow))*rho*C1*delta+lambda1*xi*exp(-xi*(B-u1-ulow))*rho*C1*b-lambda2*xi*exp(-xi*(B-u1-ulow))*rho*C2*delta+lambda2*xi*exp(-xi*(B-u1-ulow))*rho*C2*b]']);
	
	case 10
		out=mystr2sym(['[0,' ...
			'0,' ...
			'-(-kappa*B+kappa*u2+kappa*vlow-lambda1*xi*exp(-xi*vlow)*rho*C1*delta*B+lambda1*xi*exp(-xi*vlow)*rho*C1*delta*u2+lambda1*xi*exp(-xi*vlow)*rho*C1*delta*vlow+lambda1*xi*exp(-xi*vlow)*rho*C1*b*B-lambda1*xi*exp(-xi*vlow)*rho*C1*b*u2-lambda1*xi*exp(-xi*vlow)*rho*C1*b*vlow-lambda2*xi*exp(-xi*vlow)*rho*C2*delta*B+lambda2*xi*exp(-xi*vlow)*rho*C2*delta*u2+lambda2*xi*exp(-xi*vlow)*rho*C2*delta*vlow+lambda2*xi*exp(-xi*vlow)*rho*C2*b*B-lambda2*xi*exp(-xi*vlow)*rho*C2*b*u2-lambda2*xi*exp(-xi*vlow)*rho*C2*b*vlow+gamma1*B-gamma1*u2-gamma1*vlow+lambda1*eta*(B-u2-vlow)^alpha*alpha*(C1+epsilon)^beta)/(B-u2-vlow),' ...
			'-(gamma1*B-gamma1*u2-gamma1*vlow+lambda1*eta*(B-u2-vlow)^alpha*alpha*(C1+epsilon)^beta)/(B-u2-vlow)]']);
	
	case 11
		out=mystr2sym(['[(gamma1*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta)/ulow,' ...
			'(gamma2*ulow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta)/ulow,' ...
			'kappa+lambda1*xi*exp(-xi*vlow)*rho*C1*delta-lambda1*xi*exp(-xi*vlow)*rho*C1*b+lambda2*xi*exp(-xi*vlow)*rho*C2*delta-lambda2*xi*exp(-xi*vlow)*rho*C2*b,' ...
			'0]']);
	
	case 12
		out=mystr2sym(['[(gamma1*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta-ulow*kappa-ulow*lambda1*xi*exp(-xi*(B-2*ulow))*rho*C1*delta+ulow*lambda1*xi*exp(-xi*(B-2*ulow))*rho*C1*b-ulow*lambda2*xi*exp(-xi*(B-2*ulow))*rho*C2*delta+ulow*lambda2*xi*exp(-xi*(B-2*ulow))*rho*C2*b)/ulow,' ...
			'(gamma2*ulow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta-ulow*kappa-ulow*lambda1*xi*exp(-xi*(B-2*ulow))*rho*C1*delta+ulow*lambda1*xi*exp(-xi*(B-2*ulow))*rho*C1*b-ulow*lambda2*xi*exp(-xi*(B-2*ulow))*rho*C2*delta+ulow*lambda2*xi*exp(-xi*(B-2*ulow))*rho*C2*b)/ulow,' ...
			'0,' ...
			'-kappa-lambda1*xi*exp(-xi*(B-2*ulow))*rho*C1*delta+lambda1*xi*exp(-xi*(B-2*ulow))*rho*C1*b-lambda2*xi*exp(-xi*(B-2*ulow))*rho*C2*delta+lambda2*xi*exp(-xi*(B-2*ulow))*rho*C2*b]']);
	
	case 13
		out=mystr2sym(['[-(-gamma1*ulow*B+gamma1*ulow^2+gamma1*ulow*vlow-lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta*B+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta*ulow+lambda1*eta*ulow^alpha*alpha*(C1+epsilon)^beta*vlow+ulow*gamma2*B-gamma2*ulow^2-ulow*gamma2*vlow+ulow*lambda2*eta*(B-ulow-vlow)^alpha*alpha*(C2+epsilon)^beta)/ulow/(B-ulow-vlow),' ...
			'0,' ...
			'-(-kappa*B+ulow*kappa+kappa*vlow-lambda1*xi*exp(-xi*vlow)*rho*C1*delta*B+lambda1*xi*exp(-xi*vlow)*rho*C1*delta*ulow+lambda1*xi*exp(-xi*vlow)*rho*C1*delta*vlow+lambda1*xi*exp(-xi*vlow)*rho*C1*b*B-lambda1*xi*exp(-xi*vlow)*rho*C1*b*ulow-lambda1*xi*exp(-xi*vlow)*rho*C1*b*vlow-lambda2*xi*exp(-xi*vlow)*rho*C2*delta*B+lambda2*xi*exp(-xi*vlow)*rho*C2*delta*ulow+lambda2*xi*exp(-xi*vlow)*rho*C2*delta*vlow+lambda2*xi*exp(-xi*vlow)*rho*C2*b*B-lambda2*xi*exp(-xi*vlow)*rho*C2*b*ulow-lambda2*xi*exp(-xi*vlow)*rho*C2*b*vlow+gamma2*B-gamma2*ulow-gamma2*vlow+lambda2*eta*(B-ulow-vlow)^alpha*alpha*(C2+epsilon)^beta)/(B-ulow-vlow),' ...
			'-(gamma2*B-gamma2*ulow-gamma2*vlow+lambda2*eta*(B-ulow-vlow)^alpha*alpha*(C2+epsilon)^beta)/(B-ulow-vlow)]']);
	
	case 14
		out=mystr2sym(['[0,' ...
			'(ulow*gamma2*B-gamma2*ulow^2-ulow*gamma2*vlow+lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta*B-lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta*ulow-lambda2*eta*ulow^alpha*alpha*(C2+epsilon)^beta*vlow-gamma1*ulow*B+gamma1*ulow^2+gamma1*ulow*vlow-ulow*lambda1*eta*(B-ulow-vlow)^alpha*alpha*(C1+epsilon)^beta)/ulow/(B-ulow-vlow),' ...
			'-(-kappa*B+ulow*kappa+kappa*vlow-lambda1*xi*exp(-xi*vlow)*rho*C1*delta*B+lambda1*xi*exp(-xi*vlow)*rho*C1*delta*ulow+lambda1*xi*exp(-xi*vlow)*rho*C1*delta*vlow+lambda1*xi*exp(-xi*vlow)*rho*C1*b*B-lambda1*xi*exp(-xi*vlow)*rho*C1*b*ulow-lambda1*xi*exp(-xi*vlow)*rho*C1*b*vlow-lambda2*xi*exp(-xi*vlow)*rho*C2*delta*B+lambda2*xi*exp(-xi*vlow)*rho*C2*delta*ulow+lambda2*xi*exp(-xi*vlow)*rho*C2*delta*vlow+lambda2*xi*exp(-xi*vlow)*rho*C2*b*B-lambda2*xi*exp(-xi*vlow)*rho*C2*b*ulow-lambda2*xi*exp(-xi*vlow)*rho*C2*b*vlow+gamma1*B-gamma1*ulow-gamma1*vlow+lambda1*eta*(B-ulow-vlow)^alpha*alpha*(C1+epsilon)^beta)/(B-ulow-vlow),' ...
			'-(gamma1*B-gamma1*ulow-gamma1*vlow+lambda1*eta*(B-ulow-vlow)^alpha*alpha*(C1+epsilon)^beta)/(B-ulow-vlow)]']);
	
end
	
	
%-------------------------------------------------------------------------
function out=cartel2dsymSymbolicStateLagrangeMultiplier(arcid)
% returns the symbolic expressions for the Lagrangian multipliers for pure state inequality constraints for cartel2dsym
	
out=sym([]);
