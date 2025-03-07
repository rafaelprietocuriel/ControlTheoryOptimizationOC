function DXDt=cartel2dsymCanonicalSystemGeneralized(t,depvar,par,arcid)
%
% CARTEL2DSYMCANONICALSYSTEMGENERALIZED returns the canonical system together with the control dynamics
% for implicitly given controls.
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
DXDt=cartel2dsymCanonicalSystem(t,depvar,par,arcid);
DUDt=cartel2dsymImplicitControlDynamics(t,depvar,par,arcid);
	
DXDt=[DXDt;DUDt];
