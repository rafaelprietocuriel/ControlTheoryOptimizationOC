function coef=calchnf_LP(tmesh,coeff,tangent,operatoreq,phi,psi,J)
global OCMATCONT 

[t,y,z,freepar,modelpar]=OCMATCONT.drearr(tmesh,coeff);
if nargin==6
    M=calc_RHSJac(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,[]);
    M=reduceJac(M,phi,psi);
else
    M=[J psi;phi',0];
end
b=[];
b(OCMATCONT.HE.numdvariables-OCMATCONT.codimension+1,1)=1;

w=M\b;
v=M'\b;
v=v(1:OCMATCONT.HE.numdvariables-1);
w=w(1:OCMATCONT.HE.numdvariables-1);
w=w/norm(w);
v=v/(v'*w);
switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        h2=multlinear2_bvp4c(t,y,freepar,modelpar,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.odehess,OCMATCONT.bchess,w,w);
    otherwise
        h2=[];
end
if isempty(h2)
    increment=OCMATCONT.OPTIONS.increment^(3.0/4.0);
    h2=calcmultilinear2(operatoreq,w,w,tmesh,coeff,tangent,increment);
end
coef=v'*h2(:)/2.0;
%----------------------------------------------------

function J=reduceJac(J,phi,psi)
global OCMATCONT

J(:,OCMATCONT.HE.numdvariables)=[];
J(1:OCMATCONT.HE.numdvariables-OCMATCONT.codimension,OCMATCONT.HE.numdvariables-OCMATCONT.codimension+1)=psi(:); 
J(OCMATCONT.HE.numdvariables-OCMATCONT.codimension+1,1:OCMATCONT.HE.numdvariables-OCMATCONT.codimension)=phi(:).';
J(OCMATCONT.HE.numdvariables-OCMATCONT.codimension+1,OCMATCONT.HE.numdvariables-OCMATCONT.codimension+1)=0;
