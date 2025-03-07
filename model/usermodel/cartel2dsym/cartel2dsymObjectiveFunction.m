function out=cartel2dsymObjectiveFunction(t,depvar,par,arcid)
% returns the objective function for cartel2dsym
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
	
	
out=exp(-r.*t).*(-w.*s.*depvar(1,:).*depvar(2,:)-(1-w).*h.*(depvar(1,:)+depvar(2,:))-gamma1.*ctrl(1,:)-gamma2.*ctrl(2,:)-kappa.*ctrl(3,:));
