function J=cartel2dsymImplicitControlDynamicsParameterJacobian(t,depvar,par,arcid)
%
% returns the Jacobian of the implicit control dynamics with respect to the state, 
% costate and implicit control.
%
% J=CARTEL2DSYMIMPLICITCONTROLDYNAMICSPARAMETERJACOBIAN(t,depvar,par,arcid)
%
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
dim=[1,1,1,0,2,1,0,1,0,1,1,0,0,0,0]; % number of implicit controls
	
%-----------------------
% numerical calculation
%-----------------------
%numJacOpt.diffvar=3;
%numJacOpt.vectvars=[];
%if dim(arcid+1)>0
%	Jnum=numjaccsd(@cartel2dsymImplicitControlDynamics,{t,depvar,par(:),arcid},dim(arcid+1),numJacOpt);
%	Jnum=Jnum(:,1:(4+dim(arcid+1)));
%else
%	Jnum=[];
%end
%J=Jnum;
%return
%-----------------------
% analytic calculation
%-----------------------
DctrlDt=zeros(3,1);
if dim(arcid+1)
	[dudp,coord,tmp]=cartel2dsymDImplicitControlDP(t,depvar,par,arcid);
	DctrlDt(coord)=tmp;
	PUU=cartel2dsymD2LagrangianDU2(t,depvar,par,arcid);
	PUU=reshape(PUU,dim(arcid+1),dim(arcid+1));
	GP=cartel2dsymImplicitGP(t,depvar,par,arcid,DctrlDt);
	GU=cartel2dsymImplicitGU(t,depvar,par,arcid,DctrlDt);
	
	J(1:dim(arcid+1),1:21)=-PUU\(GP+GU*dudp);
else
	J=[];
end
% We used the function G=LUX*dot-X+LUU*dot-u, yielding G_X+G_u*dudx+G_{dot-u}D(dot-u(X))/DX
% with G_{dot-u}=LUU and D(dot-u(X))/DX the total derivative of dot-u with respect to X.
% The total derivative of D(dot-u(X))/DX = d(dot-u(X))/du*dU/dX+d(dot-u(X))/dX = -PUU\(G_X+G_u*dudx)
% this yields d(dot-u(X))/dX (partial derivative with respect to X) = -PUU\GX and
% d(dot-u(X))/du (partial derivative with respect to u) = -PUU\GU
