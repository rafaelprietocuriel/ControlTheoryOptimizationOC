function J=cartel2dsymCanonicalSystemGeneralizedJacobian(t,depvar,par,arcid)
%
% CARTEL2DSYMCANONICALSYSTEMGENERALIZEDGENERALIZEDJACOBIAN() returns the jacobian
% of the generalized canonical system, i.e. together with the control dynamics for
% implicitly given controls, with respect to the generalized state.
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
J=cartel2dsymCanonicalSystemJacobian(t,depvar,par,arcid);
JU=cartel2dsymImplicitControlDynamicsJacobian(t,depvar,par,arcid);
	
J=[J;JU];
