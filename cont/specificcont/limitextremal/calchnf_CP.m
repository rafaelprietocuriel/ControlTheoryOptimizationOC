function coef=calchnf_CP(tmesh,coeff,tangent,operatoreq,phi,psi,J)
global OCMATCONT 

[t,y,z,freepar,modelpar]=OCMATCONT.drearr(tmesh,coeff);
if nargin==6
    M=calc_RHSJac(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,[]);
    M=reduceJac(M,phi,psi);
else
    M=[J psi(:);phi(:)',0];
end
b=[];
b(OCMATCONT.HE.numdvariables-1,1)=1;

w=M\b;
v=M'\b;
v=v(1:OCMATCONT.HE.numdvariables-2);
w=w(1:OCMATCONT.HE.numdvariables-2);
w=w/norm(w);
v=v/(v'*w);

% see kuznetsov1999
%Y.A.~Kuznetsov, Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in Ode's},
% SIAM Journal on Numerical Analysis, 1999, 4 (36), 1104-1124
switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        h2=[];
    otherwise
        h2=[];
end
if isempty(h2)
    hessIncrement = OCMATCONT.OPTIONS.increment^(3.0/4.0);
    ten3Increment = OCMATCONT.OPTIONS.increment^(3.0/5.0);
    h2 = calcmultilinear2(operatoreq,w,w,tmesh,coeff,tangent,hessIncrement);
    f2 = v'*h2/2.0;
    h2 = [J w ; v' 0] \ [2*f2*w-h2 ; 0]; 				%-A^{INV}( B(q0,q0) - <p0,B(q0,q0)>q0 )
    h2 = h2(1:OCMATCONT.HE.numdvariablesmc);
    h3 = calcmultilinear3(operatoreq,w,w,w,tmesh,coeff,tangent,ten3Increment);
    h3 = h3 +3*calcmultilinear2(operatoreq,w,h2,tmesh,coeff,tangent,hessIncrement);	%+3B(q0,h2)
    coef = v'*h3/6.0;
end

function J=reduceJac(J,phi,psi)
global OCMATCONT

J(:,OCMATCONT.HE.numdvariables)=[];
J(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=psi(:); % remove derivative with respect to continuation parameter
J(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=phi(:).'; % remove derivative with respect to continuation parameter
J(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;
