function out=cartel2dsymEquilibriumEquation(depvar,par,arcid)
% equation for an equilibrium of the canonical system
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
r=par(1);
beta=par(2);
alpha=par(3);
rho=par(4);
theta=par(5);
omega=par(6);
delta=par(7);
h=par(8);
gamma1=par(9);
gamma2=par(10);
kappa=par(11);
eta=par(12);
ulow=par(13);
vlow=par(14);
w=par(15);
s=par(16);
b=par(17);
xi=par(18);
epsilon=par(19);
tau=par(20);
B=par(21);
	
ctrl=cartel2dsymOptimalControl(0,depvar,par,arcid);
	
out=[tau+(b+(delta-b).*exp(-xi.*ctrl(3,:))).*rho.*depvar(1,:)-eta.*ctrl(1,:).^alpha.*(depvar(1,:)+epsilon).^beta-theta.*depvar(1,:).*depvar(2,:)-omega.*depvar(1,:).^2; ...
	tau+(b+(delta-b).*exp(-xi.*ctrl(3,:))).*rho.*depvar(2,:)-eta.*ctrl(2,:).^alpha.*(depvar(2,:)+epsilon).^beta-theta.*depvar(1,:).*depvar(2,:)-omega.*depvar(2,:).^2; ...
	r.*depvar(3,:)-(-w.*s.*depvar(2,:)-(1-w).*h+depvar(3,:).*((b+(delta-b).*exp(-xi.*ctrl(3,:))).*rho-eta.*ctrl(1,:).^alpha.*(depvar(1,:)+epsilon).^beta.*beta./(depvar(1,:)+epsilon)-theta.*depvar(2,:)-2.*omega.*depvar(1,:))-depvar(4,:).*theta.*depvar(2,:)); ...
	r.*depvar(4,:)-( -w.*s.*depvar(1,:)-(1-w).*h-depvar(3,:).*theta.*depvar(1,:)+depvar(4,:).*((b+(delta-b).*exp(-xi.*ctrl(3,:))).*rho-eta.*ctrl(2,:).^alpha.*(depvar(2,:)+epsilon).^beta.*beta./(depvar(2,:)+epsilon)-theta.*depvar(1,:)-2.*omega.*depvar(2,:)))];
	
% Algebraic equations
outae=cartel2dsymAlgebraicEquation(0,depvar,par,arcid);
out=[out;outae];
