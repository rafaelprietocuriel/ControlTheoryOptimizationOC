function [lagmcc,lagmsc]=cartel2dsymLagrangeMultiplier(t,depvar,par,arcid)
% returns the Lagrangian multipliers for inequality constraints for cartel2dsym
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
lagmcc=cartel2dsymControlLagrangeMultiplier(t,depvar,par,arcid);
lagmsc=cartel2dsymStateLagrangeMultiplier(t,depvar,par,arcid);
	
%-------------------------------------------------------------------------
function out=cartel2dsymControlLagrangeMultiplier(t,depvar,par,arcid)
% returns the Lagrangian multipliers for (mixed) control inequality constraints for cartel2dsym
	
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
	
switch arcid
	case 0
		out=[repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t))];
	
	case 1
		out=[(gamma1.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./ulow; ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t))];
	
	case 2
		out=[repmat(0,1,length(t)); ...
			(gamma2.*ulow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./ulow; ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t))];
	
	case 3
		out=[repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			kappa+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b; ...
			repmat(0,1,length(t))];
	
	case 4
		out=[repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			-kappa-depvar(3,:).*xi.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*delta+depvar(3,:).*xi.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*b-depvar(4,:).*xi.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*delta+depvar(4,:).*xi.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*b];
	
	case 5
		out=[(gamma1.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./ulow; ...
			(gamma2.*ulow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./ulow; ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t))];
	
	case 6
		out=[(gamma1.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./ulow; ...
			repmat(0,1,length(t)); ...
			kappa+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b; ...
			repmat(0,1,length(t))];
	
	case 7
		out=[(gamma1.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta-ulow.*kappa-ulow.*depvar(3,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(1,:).*delta+ulow.*depvar(3,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(1,:).*b-ulow.*depvar(4,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(2,:).*delta+ulow.*depvar(4,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(2,:).*b)./ulow; ...
			repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			-kappa-depvar(3,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(1,:).*delta+depvar(3,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(1,:).*b-depvar(4,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(2,:).*delta+depvar(4,:).*xi.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(2,:).*b];
	
	case 8
		out=[repmat(0,1,length(t)); ...
			(gamma2.*ulow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./ulow; ...
			kappa+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b; ...
			repmat(0,1,length(t))];
	
	case 9
		out=[repmat(0,1,length(t)); ...
			(gamma2.*ulow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta-ulow.*kappa-ulow.*depvar(3,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(1,:).*delta+ulow.*depvar(3,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(1,:).*b-ulow.*depvar(4,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(2,:).*delta+ulow.*depvar(4,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(2,:).*b)./ulow; ...
			repmat(0,1,length(t)); ...
			-kappa-depvar(3,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(1,:).*delta+depvar(3,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(1,:).*b-depvar(4,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(2,:).*delta+depvar(4,:).*xi.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(2,:).*b];
	
	case 10
		out=[repmat(0,1,length(t)); ...
			repmat(0,1,length(t)); ...
			-(-kappa.*B+kappa.*depvar(5,:)+kappa.*vlow-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*B+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*depvar(5,:)+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*vlow+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*B-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*depvar(5,:)-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*vlow-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*B+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*depvar(5,:)+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*vlow+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*B-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*depvar(5,:)-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*vlow+gamma1.*B-gamma1.*depvar(5,:)-gamma1.*vlow+depvar(3,:).*eta.*(B-depvar(5,:)-vlow).^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./(B-depvar(5,:)-vlow); ...
			-(gamma1.*B-gamma1.*depvar(5,:)-gamma1.*vlow+depvar(3,:).*eta.*(B-depvar(5,:)-vlow).^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./(B-depvar(5,:)-vlow)];
	
	case 11
		out=[(gamma1.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./ulow; ...
			(gamma2.*ulow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./ulow; ...
			kappa+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b; ...
			repmat(0,1,length(t))];
	
	case 12
		out=[(gamma1.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta-ulow.*kappa-ulow.*depvar(3,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(1,:).*delta+ulow.*depvar(3,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(1,:).*b-ulow.*depvar(4,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(2,:).*delta+ulow.*depvar(4,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(2,:).*b)./ulow; ...
			(gamma2.*ulow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta-ulow.*kappa-ulow.*depvar(3,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(1,:).*delta+ulow.*depvar(3,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(1,:).*b-ulow.*depvar(4,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(2,:).*delta+ulow.*depvar(4,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(2,:).*b)./ulow; ...
			repmat(0,1,length(t)); ...
			-kappa-depvar(3,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(1,:).*delta+depvar(3,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(1,:).*b-depvar(4,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(2,:).*delta+depvar(4,:).*xi.*exp(-xi.*(B-2.*ulow)).*rho.*depvar(2,:).*b];
	
	case 13
		out=[-(-gamma1.*ulow.*B+gamma1.*ulow.^2+gamma1.*ulow.*vlow-depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta.*B+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta.*ulow+depvar(3,:).*eta.*ulow.^alpha.*alpha.*(depvar(1,:)+epsilon).^beta.*vlow+ulow.*gamma2.*B-gamma2.*ulow.^2-ulow.*gamma2.*vlow+ulow.*depvar(4,:).*eta.*(B-ulow-vlow).^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./ulow./(B-ulow-vlow); ...
			repmat(0,1,length(t)); ...
			-(-kappa.*B+ulow.*kappa+kappa.*vlow-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*B+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*ulow+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*vlow+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*B-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*ulow-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*vlow-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*B+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*ulow+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*vlow+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*B-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*ulow-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*vlow+gamma2.*B-gamma2.*ulow-gamma2.*vlow+depvar(4,:).*eta.*(B-ulow-vlow).^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./(B-ulow-vlow); ...
			-(gamma2.*B-gamma2.*ulow-gamma2.*vlow+depvar(4,:).*eta.*(B-ulow-vlow).^alpha.*alpha.*(depvar(2,:)+epsilon).^beta)./(B-ulow-vlow)];
	
	case 14
		out=[repmat(0,1,length(t)); ...
			(ulow.*gamma2.*B-gamma2.*ulow.^2-ulow.*gamma2.*vlow+depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta.*B-depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta.*ulow-depvar(4,:).*eta.*ulow.^alpha.*alpha.*(depvar(2,:)+epsilon).^beta.*vlow-gamma1.*ulow.*B+gamma1.*ulow.^2+gamma1.*ulow.*vlow-ulow.*depvar(3,:).*eta.*(B-ulow-vlow).^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./ulow./(B-ulow-vlow); ...
			-(-kappa.*B+ulow.*kappa+kappa.*vlow-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*B+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*ulow+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*delta.*vlow+depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*B-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*ulow-depvar(3,:).*xi.*exp(-xi.*vlow).*rho.*depvar(1,:).*b.*vlow-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*B+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*ulow+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*delta.*vlow+depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*B-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*ulow-depvar(4,:).*xi.*exp(-xi.*vlow).*rho.*depvar(2,:).*b.*vlow+gamma1.*B-gamma1.*ulow-gamma1.*vlow+depvar(3,:).*eta.*(B-ulow-vlow).^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./(B-ulow-vlow); ...
			-(gamma1.*B-gamma1.*ulow-gamma1.*vlow+depvar(3,:).*eta.*(B-ulow-vlow).^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./(B-ulow-vlow)];
	
end
	
	
%-------------------------------------------------------------------------
function out=cartel2dsymStateLagrangeMultiplier(t,depvar,par,arcid)
% returns the Lagrangian multipliers for pure state inequality constraints for cartel2dsym
	
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
	
out=[];
