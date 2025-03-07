function out=cartel2dsymUserFunction(t,depvar,par,arcid)
% the output of this file is called by userfunction(ocObj,ocElement)
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
	
%[lagmcc,lagmsc]=cartel2dsymLagrangeMultiplier(t,depvar,par,arcid);
out=[eta.*ctrl(1,:).^alpha.*(depvar(1,:)+epsilon).^beta; ...
    eta.*ctrl(2,:).^alpha.*(depvar(2,:)+epsilon).^beta];
