function out=cartel2dsymAdmissible(t,depvar,par,arcid)
% returns the values of the constraint and Lagrangian multiplier, the solution is 
% admissible if these values are non-negative.
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
constraintval=cartel2dsymConstraint(t,depvar,par,arcid);
	
lagmcc=cartel2dsymLagrangeMultiplier(t,depvar,par,arcid);
out=[constraintval;lagmcc];
