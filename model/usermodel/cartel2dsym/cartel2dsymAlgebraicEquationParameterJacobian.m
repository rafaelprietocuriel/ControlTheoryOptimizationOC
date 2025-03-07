function out=cartel2dsymAlgebraicEquationParameterJacobian(t,depvar,par,arcid)
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
		out=[0, 0, 0, -depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(1)-depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(2), 0, 0, -depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0, -1, 0, 0, 0, 0, 0, depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), -depvar(3).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(3).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(2)+depvar(4).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0];
	
	case 1
		out=[0, 0, 0, -depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(1)-depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(2), 0, 0, -depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0, -1, 0, 0, 0, 0, 0, depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), -depvar(3).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(3).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(2)+depvar(4).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0];
	
	case 2
		out=[0, 0, 0, -depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(1)-depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(2), 0, 0, -depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0, -1, 0, 0, 0, 0, 0, depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), -depvar(3).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(3).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(2)+depvar(4).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0];
	
	case 3
		out=[];
	
	case 4
		out=[0, -depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta.*log(depvar(1)+epsilon), -depvar(3).*eta.*ctrl(1).^alpha.*log(ctrl(1)).*alpha./ctrl(1).*(depvar(1)+epsilon).^beta-depvar(3).*eta.*ctrl(1).^alpha./ctrl(1).*(depvar(1)+epsilon).^beta, depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(1).*delta-depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(1).*b+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(2).*delta-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(2).*b, 0, 0, depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2), 0, -1, 0, 1, -depvar(3).*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta, 0, 0, 0, 0, -depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2), depvar(3).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta+depvar(3).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b-depvar(3).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta+depvar(4).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b-depvar(4).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b, -depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta.*beta./(depvar(1)+epsilon), 0, -depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta+depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b; ...
			0, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*log(depvar(2)+epsilon), -depvar(4).*eta.*ctrl(2).^alpha.*log(ctrl(2)).*alpha./ctrl(2).*(depvar(2)+epsilon).^beta-depvar(4).*eta.*ctrl(2).^alpha./ctrl(2).*(depvar(2)+epsilon).^beta, depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(1).*delta-depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(1).*b+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(2).*delta-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*depvar(2).*b, 0, 0, depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2), 0, 0, -1, 1, -depvar(4).*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta, 0, 0, 0, 0, -depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2), depvar(3).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta+depvar(3).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b-depvar(3).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta+depvar(4).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b-depvar(4).*xi.*(-B+ctrl(1)+ctrl(2)).*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*beta./(depvar(2)+epsilon), 0, -depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*delta+depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(1).*b-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*delta+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ctrl(2))).*rho.*depvar(2).*b];
	
	case 5
		out=[0, 0, 0, -depvar(3).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(1)-depvar(4).*(delta-b).*xi.*exp(-xi.*ctrl(3)).*depvar(2), 0, 0, -depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0, -1, 0, 0, 0, 0, 0, depvar(3).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*ctrl(3)).*rho.*depvar(2), -depvar(3).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(1)+depvar(3).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(1)-depvar(4).*(delta-b).*exp(-xi.*ctrl(3)).*rho.*depvar(2)+depvar(4).*(delta-b).*xi.*ctrl(3).*exp(-xi.*ctrl(3)).*rho.*depvar(2), 0, 0, 0];
	
	case 6
		out=[];
	
	case 7
		out=[0, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*log(depvar(2)+epsilon), -depvar(4).*eta.*ctrl(2).^alpha.*log(ctrl(2)).*alpha./ctrl(2).*(depvar(2)+epsilon).^beta-depvar(4).*eta.*ctrl(2).^alpha./ctrl(2).*(depvar(2)+epsilon).^beta, depvar(3).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*depvar(1).*delta-depvar(3).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*depvar(1).*b+depvar(4).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*depvar(2).*delta-depvar(4).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*depvar(2).*b, 0, 0, depvar(3).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2), 0, 0, -1, 1, -depvar(4).*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta, depvar(3).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*b, 0, 0, 0, -depvar(3).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2), depvar(3).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*delta+depvar(3).*xi.*(-B+ulow+ctrl(2)).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*delta-depvar(3).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*b-depvar(3).*xi.*(-B+ulow+ctrl(2)).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*b+depvar(4).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*delta+depvar(4).*xi.*(-B+ulow+ctrl(2)).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*delta-depvar(4).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*b-depvar(4).*xi.*(-B+ulow+ctrl(2)).*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*b, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*beta./(depvar(2)+epsilon), 0, -depvar(3).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*delta+depvar(3).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(1).*b-depvar(4).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*delta+depvar(4).*xi.^2.*exp(-xi.*(B-ulow-ctrl(2))).*rho.*depvar(2).*b];
	
	case 8
		out=[];
	
	case 9
		out=[0, -depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta.*log(depvar(1)+epsilon), -depvar(3).*eta.*ctrl(1).^alpha.*log(ctrl(1)).*alpha./ctrl(1).*(depvar(1)+epsilon).^beta-depvar(3).*eta.*ctrl(1).^alpha./ctrl(1).*(depvar(1)+epsilon).^beta, depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*depvar(1).*delta-depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*depvar(1).*b+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*depvar(2).*delta-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*depvar(2).*b, 0, 0, depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1)+depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2), 0, -1, 0, 1, -depvar(3).*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta, depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*delta-depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*b+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*delta-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*b, 0, 0, 0, -depvar(3).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1)-depvar(4).*xi.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2), depvar(3).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*delta+depvar(3).*xi.*(-B+ctrl(1)+ulow).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*delta-depvar(3).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*b-depvar(3).*xi.*(-B+ctrl(1)+ulow).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*b+depvar(4).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*delta+depvar(4).*xi.*(-B+ctrl(1)+ulow).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*delta-depvar(4).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*b-depvar(4).*xi.*(-B+ctrl(1)+ulow).*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*b, -depvar(3).*eta.*ctrl(1).^alpha.*alpha./ctrl(1).*(depvar(1)+epsilon).^beta.*beta./(depvar(1)+epsilon), 0, -depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*delta+depvar(3).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(1).*b-depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*delta+depvar(4).*xi.^2.*exp(-xi.*(B-ctrl(1)-ulow)).*rho.*depvar(2).*b];
	
	case 10
		out=[0, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*log(depvar(2)+epsilon)+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta.*log(depvar(1)+epsilon)./(B-ctrl(2)-vlow), -depvar(4).*eta.*ctrl(2).^alpha.*log(ctrl(2)).*alpha./ctrl(2).*(depvar(2)+epsilon).^beta-depvar(4).*eta.*ctrl(2).^alpha./ctrl(2).*(depvar(2)+epsilon).^beta+(depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*log(B-ctrl(2)-vlow).*alpha.*(depvar(1)+epsilon).^beta+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow), 0, 0, 0, 0, 0, 1, -1, 0, -depvar(4).*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta+depvar(3).*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta./(B-ctrl(2)-vlow), 0, (-gamma1-depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.^2./(B-ctrl(2)-vlow).*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow)+(gamma1.*B-gamma1.*ctrl(2)-gamma1.*vlow+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow).^2, 0, 0, 0, 0, -depvar(4).*eta.*ctrl(2).^alpha.*alpha./ctrl(2).*(depvar(2)+epsilon).^beta.*beta./(depvar(2)+epsilon)+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta.*beta./(depvar(1)+epsilon)./(B-ctrl(2)-vlow), 0, (gamma1+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.^2./(B-ctrl(2)-vlow).*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow)-(gamma1.*B-gamma1.*ctrl(2)-gamma1.*vlow+depvar(3).*eta.*(B-ctrl(2)-vlow).^alpha.*alpha.*(depvar(1)+epsilon).^beta)./(B-ctrl(2)-vlow).^2];
	
	case 11
		out=[];
	
	case 12
		out=[];
	
	case 13
		out=[];
	
	case 14
		out=[];
	
end
