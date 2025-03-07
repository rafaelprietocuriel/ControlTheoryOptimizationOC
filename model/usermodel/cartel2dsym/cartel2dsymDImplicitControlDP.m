function [dudp,coord,dudt]=cartel2dsymDImplicitControlDP(t,depvar,par,arcid)
%
% CARTEL2DSYMDIMPLICITCONTROLDP returns the derivative of the generalized control with respect to the generalized state.
%
% DUDX=CARTEL2DSYMDIMPLICITCONTROLDX(t,depvar,par,arcid)
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
[luu,coord]=cartel2dsymD2LagrangianDU2(t,depvar,par,arcid);
dim=[1,1,1,0,2,1,0,1,0,1,1,0,0,0,0]; % number of implicit controls
if dim(arcid+1)
	lux=cartel2dsymD2LagrangianDUDX(t,depvar,par,arcid);
	lup=cartel2dsymD2LagrangianDUDP(t,depvar,par,arcid);
	dxdt=cartel2dsymCanonicalSystem(t,depvar,par,arcid);
	
	luu=reshape(luu,dim(arcid+1),dim(arcid+1),[]);
	lux=reshape(lux,dim(arcid+1),4,[]);
	lup=reshape(lup,dim(arcid+1),21,[]);
	dudp=zeros(dim(arcid+1),21,length(t));
	dudt=zeros(dim(arcid+1),length(t));
	for ii=1:length(t)
		dudp(:,:,ii)=-luu(:,:,ii)\lup(:,:,ii);
		dudt(:,ii)=(-luu(:,:,ii)\lux(:,:,ii))*dxdt(1:4,ii);
	end
else
	dudp=[];
	dudt=[];
end
% We used the implicit derivative. From LU=0 we find LUX+LUU*dudx=0 and LUX*dot-X+LUU*dot-u=0 yielding
% dot-u=-LUU^{-1}*LUX*dot-X.
