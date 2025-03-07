function [dudt,dudx,coord]=cartel2dsymImplicitControlDynamics(t,depvar,par,arcid)
%
% CARTEL2DSYMIMPLICITCONTROLDYNAMICS returns the dynamics of the implicit control and the deriviative with respect to the generalized state.
%
% DUDT=CARTEL2DSYMIMPLICITCONTROLDYNAMICS(t,depvar,par,arcid)
%
% [DUDT,DUDX,COORD]=CARTEL2DSYMIMPLICITCONTROLDYNAMICS(t,depvar,par,arcid)
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
[luu,coord]=cartel2dsymD2LagrangianDU2(t,depvar,par,arcid);
dim=[1,1,1,0,2,1,0,1,0,1,1,0,0,0,0]; % number of implicit controls
if dim(arcid+1)
	lux=cartel2dsymD2LagrangianDUDX(t,depvar,par,arcid);
	dxdt=cartel2dsymCanonicalSystem(t,depvar,par,arcid);
	
	luu=reshape(luu,dim(arcid+1),dim(arcid+1),[]);
	lux=reshape(lux,dim(arcid+1),4,[]);
	dudx=zeros(dim(arcid+1),4,length(t));
	dudt=zeros(dim(arcid+1),length(t));
	for ii=1:length(t)
		dudx(:,:,ii)=-luu(:,:,ii)\lux(:,:,ii);
		dudt(:,ii)=dudx(:,:,ii)*dxdt(1:4,ii);
	end
else
	dudt=[];
	dudx=[];
end
%dudt2=cartel2dsymExplicitControlDynamics(t,depvar,par,arcid);
