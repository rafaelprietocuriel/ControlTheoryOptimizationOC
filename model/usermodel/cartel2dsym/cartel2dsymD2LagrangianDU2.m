function [out,coord]=cartel2dsymD2LagrangianDU2(t,depvar,par,arcid)
%
% CARTEL2DSYMD2LAGRANGIANDU2 returns the Hessian of the Lagrangefunction with respect to the implicit
% controls.
%
% LU2=CARTEL2DSYMD2LAGRANGIANDU2(t,depvar,par,arcid) let m=#implicit controls
% n=#time grid, then LU2 is of size m*m x n.
%
% [LU2,COORD]=CARTEL2DSYMD2LAGRANGIANDU2(t,depvar,par,arcid) argument COORD returns the
% index of the implicit control.
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
	
switch arcid
	case 0
		out=depvar(3,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(1,:)+depvar(4,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(2,:);
		coord=3;
	
	case 1
		out=depvar(3,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(1,:)+depvar(4,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(2,:);
		coord=3;
	
	case 2
		out=depvar(3,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(1,:)+depvar(4,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(2,:);
		coord=3;
	
	case 3
		out=[];
		coord=[];
	
	case 4
		out=[-depvar(3,:).*eta.*depvar(5,:).^alpha.*alpha.^2./depvar(5,:).^2.*(depvar(1,:)+epsilon).^beta+depvar(3,:).*eta.*depvar(5,:).^alpha.*alpha./depvar(5,:).^2.*(depvar(1,:)+epsilon).^beta+depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*b+depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*b; ...
			depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*b+depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*b; ...
			depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*b+depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*b; ...
			-depvar(4,:).*eta.*depvar(6,:).^alpha.*alpha.^2./depvar(6,:).^2.*(depvar(2,:)+epsilon).^beta+depvar(4,:).*eta.*depvar(6,:).^alpha.*alpha./depvar(6,:).^2.*(depvar(2,:)+epsilon).^beta+depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(1,:).*b+depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-depvar(6,:))).*rho.*depvar(2,:).*b];
		coord=[1,2];
	
	case 5
		out=depvar(3,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(1,:)+depvar(4,:).*(delta-b).*xi.^2.*exp(-xi.*depvar(5,:)).*rho.*depvar(2,:);
		coord=3;
	
	case 6
		out=[];
		coord=[];
	
	case 7
		out=-depvar(4,:).*eta.*depvar(5,:).^alpha.*alpha.^2./depvar(5,:).^2.*(depvar(2,:)+epsilon).^beta+depvar(4,:).*eta.*depvar(5,:).^alpha.*alpha./depvar(5,:).^2.*(depvar(2,:)+epsilon).^beta+depvar(3,:).*xi.^2.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.^2.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(1,:).*b+depvar(4,:).*xi.^2.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.^2.*exp(-xi.*(B-ulow-depvar(5,:))).*rho.*depvar(2,:).*b;
		coord=2;
	
	case 8
		out=[];
		coord=[];
	
	case 9
		out=-depvar(3,:).*eta.*depvar(5,:).^alpha.*alpha.^2./depvar(5,:).^2.*(depvar(1,:)+epsilon).^beta+depvar(3,:).*eta.*depvar(5,:).^alpha.*alpha./depvar(5,:).^2.*(depvar(1,:)+epsilon).^beta+depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(1,:).*delta-depvar(3,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(1,:).*b+depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(2,:).*delta-depvar(4,:).*xi.^2.*exp(-xi.*(B-depvar(5,:)-ulow)).*rho.*depvar(2,:).*b;
		coord=1;
	
	case 10
		out=-depvar(4,:).*eta.*depvar(5,:).^alpha.*alpha.^2./depvar(5,:).^2.*(depvar(2,:)+epsilon).^beta+depvar(4,:).*eta.*depvar(5,:).^alpha.*alpha./depvar(5,:).^2.*(depvar(2,:)+epsilon).^beta+(-gamma1-depvar(3,:).*eta.*(B-depvar(5,:)-vlow).^alpha.*alpha.^2./(B-depvar(5,:)-vlow).*(depvar(1,:)+epsilon).^beta)./(B-depvar(5,:)-vlow)+(gamma1.*B-gamma1.*depvar(5,:)-gamma1.*vlow+depvar(3,:).*eta.*(B-depvar(5,:)-vlow).^alpha.*alpha.*(depvar(1,:)+epsilon).^beta)./(B-depvar(5,:)-vlow).^2;
		coord=2;
	
	case 11
		out=[];
		coord=[];
	
	case 12
		out=[];
		coord=[];
	
	case 13
		out=[];
		coord=[];
	
	case 14
		out=[];
		coord=[];
	
end
