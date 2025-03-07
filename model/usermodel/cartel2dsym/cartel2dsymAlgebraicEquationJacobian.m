function out=cartel2dsymAlgebraicEquationJacobian(t,depvar,par,arcid)
% returns the Jacobian of algebraic equations for the implicitly given controls
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
	
switch arcid
	case 0
		out=[-depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1), -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), depvar(3).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(2)];
	
	case 1
		out=[-depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1), -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), depvar(3).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(2)];
	
	case 2
		out=[-depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1), -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), depvar(3).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(2)];
	
	case 3
		out=[];
	
	case 4
		out=[-depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta.*beta./(depvar(1)+epsilon)+depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*delta-depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*b, depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*delta-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*b, -eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta+xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b, xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b, -depvar(3).*eta.*ctrl(1).^alpha.*alpha.^2./ctrl(1).^2.*(depvar(1)+epsilon).^beta+depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).^2.*(depvar(1)+epsilon).^beta+depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b, depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b; ...
			depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*delta-depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*b, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*beta./(depvar(2)+epsilon)+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*delta-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*b, xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b, -eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta+xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b, depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b, -depvar(4).*eta.*ctrl(2).^alpha.*alpha.^2./ctrl(2).^2.*(depvar(2)+epsilon).^beta+depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).^2.*(depvar(2)+epsilon).^beta+depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b];
	
	case 5
		out=[-depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho, -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1), -(delta-b).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), depvar(3).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*(delta-b).*xi.^2.*exp(-xi.*ctrl(3)).*rho.*depvar(2)];
	
	case 6
		out=[];
	
	case 7
		out=[depvar(3).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*delta-depvar(3).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*b, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*beta./(depvar(2)+epsilon)+depvar(4).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*delta-depvar(4).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*b, xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*delta-xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*b, -eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta+xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*delta-xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*b, -depvar(4).*eta.*ctrl(2).^alpha.*alpha.^2./ctrl(2).^2.*(depvar(2)+epsilon).^beta+depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).^2.*(depvar(2)+epsilon).^beta+depvar(3).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*b];
	
	case 8
		out=[];
	
	case 9
		out=[-depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta.*beta./(depvar(1)+epsilon)+depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*delta-depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*b, depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*delta-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*b, -eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta+xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*delta-xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*b, xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*delta-xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*b, -depvar(3).*eta.*ctrl(1).^alpha.*alpha.^2./ctrl(1).^2.*(depvar(1)+epsilon).^beta+depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).^2.*(depvar(1)+epsilon).^beta+depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*b];
	
	case 10
		out=[depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta.*beta./(depvar(1)+epsilon)./(B-ctrl(2)-vlow), -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*beta./(depvar(2)+epsilon), eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta./(B-ctrl(2)-vlow), -eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta, -depvar(4).*eta.*ctrl(2).^alpha.*alpha.^2./ctrl(2).^2.*(depvar(2)+epsilon).^beta+depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).^2.*(depvar(2)+epsilon).^beta+(-gamma1-depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.^2./(B-ctrl(2)-vlow).*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow)+(gamma1.*B-gamma1.*ctrl(2)-gamma1.*vlow+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow).^2];
	
	case 11
		out=[];
	
	case 12
		out=[];
	
	case 13
		out=[];
	
	case 14
		out=[];
	
end
