function out=cartel2dsymD2HamiltonianDX2(t,depvar,par,arcid)
% returns the Hessian of the Hamiltonian with respect to the state and costate variable(s).
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
	
out=[depvar(3)*(-eta*ctrl(1)^alpha*(depvar(1)+epsilon)^beta*beta^2/(depvar(1)+epsilon)^2+eta*ctrl(1)^alpha*(depvar(1)+epsilon)^beta*beta/(depvar(1)+epsilon)^2-2*omega), -w*s-depvar(3)*theta-depvar(4)*theta, (b+(delta-b)*exp(-xi*ctrl(3)))*rho-eta*ctrl(1)^alpha*(depvar(1)+epsilon)^beta*beta/(depvar(1)+epsilon)-theta*depvar(2)-2*omega*depvar(1), -theta*depvar(2); ...
	-w*s-depvar(3)*theta-depvar(4)*theta, depvar(4)*(-eta*ctrl(2)^alpha*(depvar(2)+epsilon)^beta*beta^2/(depvar(2)+epsilon)^2+eta*ctrl(2)^alpha*(depvar(2)+epsilon)^beta*beta/(depvar(2)+epsilon)^2-2*omega), -theta*depvar(1), (b+(delta-b)*exp(-xi*ctrl(3)))*rho-eta*ctrl(2)^alpha*(depvar(2)+epsilon)^beta*beta/(depvar(2)+epsilon)-theta*depvar(1)-2*omega*depvar(2); ...
	(b+(delta-b)*exp(-xi*ctrl(3)))*rho-eta*ctrl(1)^alpha*(depvar(1)+epsilon)^beta*beta/(depvar(1)+epsilon)-theta*depvar(2)-2*omega*depvar(1), -theta*depvar(1), 0, 0; ...
	-theta*depvar(2), (b+(delta-b)*exp(-xi*ctrl(3)))*rho-eta*ctrl(2)^alpha*(depvar(2)+epsilon)^beta*beta/(depvar(2)+epsilon)-theta*depvar(1)-2*omega*depvar(2), 0, 0];
