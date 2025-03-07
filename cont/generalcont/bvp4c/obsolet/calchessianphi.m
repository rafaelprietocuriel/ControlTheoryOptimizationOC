function dgdY=calchessianphi(t,y,freepar,modelpar,ode,bc,odejac,bcjac,odehess,bchess,reigvec,leigvec)

global OCMATCONT OCBVP
increment=1e-5;
counter=0;
  %numJacOpt.diffvar=2;
for ii=1:OCMATCONT.HE.numdvariablesmc
    counter=counter+1;
    y1 = y(:); y1(ii) = y1(ii)-increment;
    y2 = y(:); y2(ii) = y2(ii)+increment;
    J2=calc_RHSJac(t,y2(OCMATCONT.HE.DDATA.meshvalcoord),[],freepar,modelpar,ode,bc,[],odejac,bcjac,[]);
    J1=calc_RHSJac(t,y1(OCMATCONT.HE.DDATA.meshvalcoord),[],freepar,modelpar,ode,bc,[],odejac,bcjac,[]);
    J1(:,OCMATCONT.HE.numdvariables+(-1:0))=[]; % remove derivative with respect to continuation parameter
    J1(OCMATCONT.HE.numdvariables-1,:)=[]; % remove derivative with respect to continuation parameter
    J2(:,OCMATCONT.HE.numdvariables+(-1:0))=[]; % remove derivative with respect to continuation parameter
    J2(OCMATCONT.HE.numdvariables-1,:)=[]; % remove derivative with respect to continuation parameter
    dgdY(ii) = -reigvec'*((J2)-(J1))*leigvec/(2*increment);
    %numJac=mynumjaccsd(@jacobianvect,{t,y,[],freepar,modelpar,ode,bc,[],odejac,bcjac,[]},numel(y)^2,numJacOpt,ii);
end
for ii=1:2
    counter=counter+1;
    freepar1 = freepar; freepar1(ii) = freepar1(ii)-increment;
    freepar2 = freepar; freepar2(ii) = freepar2(ii)+increment;
    J2=calc_RHSJac(t,y,[],freepar2,modelpar,ode,bc,[],odejac,bcjac,[]);
    J1=calc_RHSJac(t,y,[],freepar1,modelpar,ode,bc,[],odejac,bcjac,[]);
    J1(:,OCMATCONT.HE.numdvariables+(-1:0))=[]; % remove derivative with respect to continuation parameter
    J1(OCMATCONT.HE.numdvariables-1,:)=[]; % remove derivative with respect to continuation parameter
    J2(:,OCMATCONT.HE.numdvariables+(-1:0))=[]; % remove derivative with respect to continuation parameter
    J2(OCMATCONT.HE.numdvariables-1,:)=[]; % remove derivative with respect to continuation parameter
    dgdY(counter) = -reigvec'*((J2)-(J1))*leigvec/(2*increment);
end

%  numJacOpt.vectvars=[];
%  out2=reshape(numJac,4,4,4);
% 
function J=jacobianvect(t,y,z,freepar,modelpar,ode,bc,ic,odejac,bcjac,icjac)
global OCMATCONT OCBVP

J=shallowlake2DCanonicalSystemJacobian(t,depvar,par,arcid);
J=J(:);

%-----------------------------------------------------------------
function [dFdy,nfevals,nfcalls]=mynumjaccsd(fun,Fargs,nF,options,ii)
% NUMJACCSD    Complex Step Jacobian
% based on 'jacobiancsd' by Yi Cao at Cranfield University, 02/01/2008 and
% uses the structure of odenumjac.
%
% J = NUMJACCSD(F,FARGS,NF,OPTIONS) returns the numerical (NF x N) Jacobian
% matrix of a NF-vector function, F(FARGS{:}) at the reference point, 
% Y=FARGS{OPTIONS.DIFFVAR} (N-vector). 
%
% The structure OPTIONS must have the following fields: DIFFVAR, VECTVARS.
% The field OPTIONS.DIFFVAR is the index of the  differentiation variable,
% Y = FARGS{DIFFVAR}. For a function F(t,x), set DIFFVAR to 1 to compute
% DF/Dt, or to 2 to compute DF/Dx. Set OPTIONS.VECTVAR to the indices of
% vectorized arguments: VECTVAR = [2] indicates that  F(t,[x1 y2 ...])
% returns [F(t,x1) F(t,x2) ...], while VECTVAR = [1,2] indicates that F([t1
% t2 ...],[x1 x2 ...]) returns [F(t1,x1) F(t2,x2) ...].
%   
% [DFDY,NFEVALS,NFCALLS] = ODENUMJAC(...) returns the number of values
% F(FARGS{:}) computed while forming dFdY (NFEVALS) and the number of calls
% to the function F (NFCALLS). If F is not vectorized, NFCALLS equals
% NFEVALS.  
% Dieter Grass

% Options
diffvar = options.diffvar; 

% The differentiation variable.
y  = Fargs{diffvar};
classY = class(y);

ny = 1;

dFdy=zeros(nF,ny);                   % allocate memory for the Jacobian matrix
del = (y + ny*eps(classY)*i) - y;
h=imag(del);
ydel = y(:,ones(1,ny)) + diag(del);
dFdy(:,ii)=imag(fun(Fargs{1:diffvar-1},ydel(:,ii),Fargs{diffvar+1:end}))/h(ii);     % complex step differentiation
nfcalls = ny;                       % stats
nfevals = ny;                         % stats (at least one per loop)
