function out=cartel2dsymD2HamiltonianDu2(t,depvar,par,arcid)
% returns the Hessian of the Hamiltonian with respect to the control variable(s).
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
	
ctrl=cartel2dsymOptimalControl(t,depvar,par,arcid);
	
[lagmcc,lagmsc]=cartel2dsymLagrangeMultiplier(t,depvar,par,arcid);
	
out=[-depvar(3)*eta*ctrl(1)^alpha*alpha^2/ctrl(1)^2*(depvar(1)+epsilon)^beta+depvar(3)*eta*ctrl(1)^alpha*alpha/ctrl(1)^2*(depvar(1)+epsilon)^beta, 0, 0; ...
	0, -depvar(4)*eta*ctrl(2)^alpha*alpha^2/ctrl(2)^2*(depvar(2)+epsilon)^beta+depvar(4)*eta*ctrl(2)^alpha*alpha/ctrl(2)^2*(depvar(2)+epsilon)^beta, 0; ...
	0, 0, depvar(3)*(delta-b)*xi^2*exp(-xi*ctrl(3))*rho*depvar(1)+depvar(4)*(delta-b)*xi^2*exp(-xi*ctrl(3))*rho*depvar(2)];
