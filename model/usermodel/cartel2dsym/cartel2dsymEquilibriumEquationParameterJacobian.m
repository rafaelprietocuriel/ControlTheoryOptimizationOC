function out=cartel2dsymEquilibriumEquationParameterJacobian(depvar,par,arcid)
% returns the Jacobian of the equilibrium equations
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
	
out=cartel2dsymCanonicalSystemParameterJacobian(0,depvar,par,arcid);
	
% Jacobian of algebraic equations
outae=cartel2dsymAlgebraicEquationParameterJacobian(0,depvar,par,arcid);
out=[out;outae];
