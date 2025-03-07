function out=limitextremal2lc()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{28}{1}=@ode;
out{28}{2}=@odejacobian;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;
out{5}{3}=@odehess;
out{5}{4}=@bchess;
out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
out{11}=@operatorpfrechet;
out{12}{1}=@residual;
out{12}{2}=@maxresidual;
out{12}{3}=@verifyresidual;
out{12}{4}=@prepare4meshadapt;
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{17}=@adapt;
out{18}=@meshadaptation;
out{20}=@workspaceadapt;
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
out{25}=@deval;
out{26}=@saveintermediate;
out{27}=@datapath;
out{30}=@printcontinuation;
out{31}=@test;

function res=operatoreq(tmesh,coeff,tangent,odefun,bcfun,icfunc)
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
res=calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
% append g (minimally augmented system, e.g. Kuzentsov 1998, p. 502f)
M=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfunc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
M(:,OCMATCONT.HE.numdvariables)=[];
if OCMATCONT.HE.numdvariables-2~=length(OCMATLSC.LAE_phi)
    b = []; b(OCMATCONT.HE.numdvariables-1,1)=1;
    phi = M'\b;
    phi=phi(1:OCMATCONT.HE.numdvariables-2)';
    phi=phi/norm(phi);
    psi = M\b;
    psi=psi(1:OCMATCONT.HE.numdvariables-2)';
    psi=psi/norm(psi);

    OCMATLSC.LAE_phi=phi(:);
    OCMATLSC.LAE_psi=psi(:);
    OCMATLSC.LAE_new_phi=phi(:);
    OCMATLSC.LAE_new_psi=psi(:);
end

M(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=OCMATLSC.LAE_phi'; % append right eigenvector
M(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=OCMATLSC.LAE_psi; % append left eigenvector
M(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;

%M=[J  phi]
%  [psi' 0]
% the function g is defined as
%       M*[w;g]=[0;1], w\in\R^N, (*)
% subsequently we set s=[w;g] and b=[0;1].
% If s is a solution of eq. (*) then J*w+g*phi=0. Thus, for g=0, J*w=0 and w is
% a right eigenvector of J for the eigenvalue zero.
% If s is a solution of the transposed eq. (*) then v*J'+g*psi=0, with
% v=w'. Thus, for g=0 v is a left eigenvector of J for the eigenvalue zero.

% update left and right eigenvector alternating to keep M non-singular
% The previous consideration shows that these are the according
% eigenvectors for a solution g=0.
b = []; b(OCMATCONT.HE.numdvariables-1,1)=1;
if OCMATLSC.LAE_switch==1
    % left eigenvector
    s = M'\b;
    OCMATLSC.LAE_new_psi = s(1:OCMATCONT.HE.numdvariables-2)';
else
    % right eigenvector
    s = M\b;
    OCMATLSC.LAE_new_phi= s(1:OCMATCONT.HE.numdvariables-2)';
end
res(OCMATCONT.HE.numdvariables-1,1)=s(OCMATCONT.HE.numdvariables-1);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
M=J;
M(:,OCMATCONT.HE.numdvariables)=[];
M(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=OCMATLSC.LAE_phi'; % remove derivative with respect to continuation parameter
M(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=OCMATLSC.LAE_psi; % remove derivative with respect to continuation parameter
M(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;
% calculate g'
b=[]; 
b(OCMATCONT.HE.numdvariables-1,1)=1;
w=M\b;
v=M'\b;
dgdy=calchessianphi_bvp4c(tmesh,y,freepar,modelpar,odefun,bcfun,odejacfun,bcjacfun,@odehessian,@bchessian,v(1:OCMATCONT.HE.numdvariables-2),w(1:OCMATCONT.HE.numdvariables-2));
J(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables)=dgdy;

function opt=options(opt)

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%----------------------------------------------------------------
function [res,varargout]=residual(tmesh,coeff,rhs,odefun)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        [res resindex]=residual_sbvpoc(t,y,z,freepar,modelpar,@ode,rhs);
        varargout{1}=resindex;
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
    case 'bvp4c'
        res=residual_bvp4c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
end

%----------------------------------------------------------------
function res=maxresidual(tmesh,coeff)
global OCMATCONT

switch OCMATCONT.bvpmethod
    case 'bvp5c'
        [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
        res=maxresidual_bvp5c(t,y);
    otherwise
        res=0;
end

%----------------------------------------------------------------
function b=verifyresidual(maxres)
global OCMATCONT

b=1;
if OCMATCONT.OPTIONS.meshadaptation
    switch OCMATCONT.bvpmethod
        case 'bvp5c'
            b=verifyresidual_bvp5c(maxres);
        otherwise
            b=maxres < OCMATCONT.OPTIONS.meshadaptreltol;
    end
end

%----------------------------------------------------------------
function flag=prepare4meshadapt(tmesh,coeff)
global OCMATCONT

switch OCMATCONT.bvpmethod
    case 'bvp5c'
        [t,y]=drearr(tmesh,coeff);
        prepare4meshadapt_bvp5c(t,y);
        flag=1;
    otherwise
        flag=1;
end

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    dxdt=zeros(size(depvar));
    OCBVP.multiarccalc=1;
    for ii=1:OCMATCONT.HE.numarc
        depvarint=depvar(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:);
        dxdt(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:)=ode(s,depvarint,ii,freepar,modelpar);
    end
    OCBVP.multiarccalc=0;
    return
end
arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.truncationtimecoord)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATLSC.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATLSC.objectivevaluecalc
    dxdt(OCMATLSC.objectivevaluecoord,:)=dtds*OCMATLSC.objectivefunction(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    J=zeros(OCBVP.n);
    Jpar=zeros(OCBVP.n,OCBVP.npar);
    OCBVP.multiarccalc=1;
    for ii=1:OCMATCONT.HE.numarc
        idx=OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii);
        depvarint=depvar(idx,:);
        [J(idx,idx) Jpar(idx,1:OCBVP.npar)]=odejac(s,depvarint,ii,freepar,modelpar);
    end
    OCBVP.multiarccalc=0;
    return
end
arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.truncationtimecoord)];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
if OCMATLSC.inftimetransformation && arc==OCMATCONT.HE.numarc
    t=arctime(arc)+log(arc-s)/OCMATLSC.inftimetransformation;
    dtds=-1./(OCMATLSC.inftimetransformation*(arc-s));
else
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    dtds=diffarctime(arc);
end
J=OCMATLSC.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATLSC.objectivevaluecalc
    J=[J; ...
        dtds*OCMATLSC.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
dxdt=OCMATLSC.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATCONT.HE.numarc>1
    if ~OCMATLSC.autonomous
        Jt=OCMATLSC.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATLSC.objectivevaluecalc
        dxdt(OCMATLSC.objectivevaluecoord,:)=OCMATLSC.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATLSC.objectivevaluecoord,:)=OCMATLSC.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLSC.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLSC.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        if OCMATLSC.inftimetransformation
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=0;
            %Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=-1./((exp(OCMATLSC.inftimetransformation*arctime(arc))*(arc-s)))*dxdt;
        else
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLSC.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    end
end
Jpar(1:OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATLSC.truncationtimecoord)=dxdt;

%-------------------------------------------------------------------------
function [Htot Hpar]=odehessian(s,dynVar,arc,freepar,modelpar)
global OCMATLSC OCMATCONT
arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.truncationtimecoord)];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
if OCMATLSC.inftimetransformation & arc==OCMATCONT.HE.numarc
    t=arctime(arc)+log(arc-s)/OCMATLSC.inftimetransformation;
    dtds=-1./(OCMATLSC.inftimetransformation*(arc-s));
else
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    dtds=diffarctime(arc);
end
H=OCMATLSC.canonicalsystemhessian(t,dynVar,modelpar,arcarg);
H=dtds*H;
Htot=zeros(OCMATCONT.DOMAINDDATA(arc).numode,OCMATCONT.HE.numparameter-2+OCMATCONT.DOMAINDDATA(arc).numode,OCMATCONT.DOMAINDDATA(arc).numode);
Htot(1:OCMATCONT.DOMAINDDATA(arc).numode,1:OCMATCONT.DOMAINDDATA(arc).numode,1:OCMATCONT.DOMAINDDATA(arc).numode)=H;
Hpar=zeros(OCMATCONT.DOMAINDDATA(arc).numode,OCMATCONT.HE.numparameter-2+OCMATCONT.DOMAINDDATA(arc).numode,OCMATCONT.HE.numparameter);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    res=bc(depvara,depvarb,freepar,modelpar);
    OCBVP.multiarccalc=0;
    return
end
switchtimes=freepar(OCMATLSC.switchtimecoord);
if OCMATLSC.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATLSC.jumpcostateindex)=freepar(OCMATLSC.jumpcostatecoord);
end
resconnec=[];

resinit=depvara([OCMATLSC.freecoordinate OCMATLSC.contcoordinate],1)-freepar([OCMATLSC.freecoordinateindex end]);
resasym=OCMATLSC.bcasymptotic(depvarb,OCMATLSC.asymptoticmatrix,OCMATLSC.saddlepoint);
if OCMATLSC.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end

if OCMATLSC.stateconstraint
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATLSC.reset(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATLSC.guard(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
else
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATLSC.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATLSC.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
end
res=[resinit;resconnec;resasym];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjacobian(depvara,depvarb,freepar,modelpar)
global OCMATLSC OCMATCONT
Ja=[];
Jb=[];
Jpar=[];


%-------------------------------------------------------------------------
function [Ha Hb Hpar]=bchessian(depvara,depvarb,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
Ha=zeros(OCBVP.nBCs,2*OCBVP.numode+OCBVP.nparmcod,OCBVP.numode);
Hb=Ha;
Hpar=zeros(OCBVP.nBCs,2*OCBVP.numode+OCBVP.nparmcod,OCBVP.npar);

function [tmeshnew,coeffnew,tangentnew]=meshadaptation(tmesh,coeff,tangent,res,canRemovePoints,ode)
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp4c'
        OCMATCONT.HE.DDATA.meshvalcoordold=OCMATCONT.HE.DDATA.meshvalcoord;
        OCMATCONT.HE.parametermcodcoordold=OCMATCONT.HE.parametermcodcoord;
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp4c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
        if length(coeffnew)~=length(coeff) || any(tmesh-tmeshnew)
            [t,y,z,freepar,modelpar]=drearr(tmeshnew,coeffnew);
            phi=interp1(tmesh,OCMATLSC.LAE_phi(OCMATCONT.HE.DDATA.meshvalcoordold)',tmeshnew).';
            psi=interp1(tmesh,OCMATLSC.LAE_psi(OCMATCONT.HE.DDATA.meshvalcoordold)',tmeshnew).';
            phi=phi/norm(phi);
            psi=psi/norm(psi);
            calc_RHS(t,y,z,freepar,modelpar,ode,@bc,[]);
            %append g (minimally augmented system, e.g. Kuzentsov 1998, p. 502f)
            M=calc_RHSJac(t,y,z,freepar,modelpar,ode,@bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
            M(:,OCMATCONT.HE.numdvariables)=[];
            M(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=[phi(:)' OCMATLSC.LAE_phi(OCMATCONT.HE.parametermcodcoordold)]; % append right eigenvector
            M(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=[psi(:)' OCMATLSC.LAE_psi(OCMATCONT.HE.parametermcodcoordold)]; % append left eigenvector
            M(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;

            b = []; b(OCMATCONT.HE.numdvariables-1,1)=1;
            phi = M'\b;
            phi=phi(1:OCMATCONT.HE.numdvariables-2)';
            phi=phi/norm(phi);
            psi = M\b;
            psi=psi(1:OCMATCONT.HE.numdvariables-2)';
            psi=psi/norm(psi);

            OCMATLSC.LAE_phi=phi(:);
            OCMATLSC.LAE_psi=psi(:);
            OCMATLSC.LAE_new_phi=phi(:);
            OCMATLSC.LAE_new_psi=psi(:);
        end
end


%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATLSC OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
    %arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1)); 
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATLSC.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATLSC.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
    end
    violationmat=constr<-OCMATCONT.OPTIONS.admissibletol;
    if any(violationmat(:))
        counter=counter+1;
        [rows cols]=find(violationmat);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows=rows;
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=constr;
        infoS(counter).minval=min(constr(:));
        b=min([b infoS(counter).minval]);
    end
    % test arctimes
    violationmat=diffarctime(arc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    if any(violationmat)
        counter=counter+1;
        cols=find(violationmat);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows='switchtime';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=diffarctime(arc);
        infoS(counter).minval=min(diffarctime(:));
        b=min([b infoS(counter).minval]);
    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATLSC OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
figure(1)
h=OCMATLSC.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)


function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT
idx=[];
if isempty(coeff)
    return
end
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATLSC
if nargin > 3
    s=varargin{4};
    s.data.V=varargin{3};%AsymVar.V;
    varargout{3}=s;
end

varargout{2}=nan;
% all done succesfully
varargout{1}=0;
%-------------------------------------------------------------
function [out, failed]=testfunc(id,tmesh,coeff,tangent)
global OCMATCONT OCMATLSC

out(1)=0;
failed=[];
b=[];
b(OCMATCONT.HE.numdvariables-1,1)=1;

for ii=id
    lastwarn('');
    switch ii
        case 1 % ECP
            switch OCMATCONT.bvpmethod
                case 'sbvpoc'
                case 'bvp6c'
                case 'bvp4c'
                    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                    M=calc_RHSJac(t,y,z,freepar,modelpar,@ode,@bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,[]);
                    M(:,OCMATCONT.HE.numdvariables)=[];
                    M(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=OCMATLSC.LAE_phi'; % remove derivative with respect to continuation parameter
                    M(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=OCMATLSC.LAE_psi; % remove derivative with respect to continuation parameter
                    M(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;
                    w=M\b;
                    v=M'\b;
                    out(1)=calchnf_bvp4c(tmesh,y,freepar,modelpar,@ode,@bc,OCMATCONT.odejac,OCMATCONT.bcjac,@odehessian,@bchessian,v(1:OCMATCONT.HE.numdvariables-2),w(1:OCMATCONT.HE.numdvariables-2));
                    %out(1)=tangent(end);
                otherwise
                    error('No such testfunction');
            end
    end
    if 0%~isempty(lastwarn)
        failed=[failed ii];
    end

end
%----------------------------------------------------------------
function [out, failed]=targetvaluefunc(id,tmesh,coeff,tangent)
global OCMATCONT OCMATLSC

failed=[];
for ii=id
    switch ii
        case 1
            out=OCMATLSC.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATLSC OCBVP
dataadaptation(tmesh);
[t,y,z,freepar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATLSC.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.truncationtimecoord)];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATLSC.initialtime;
out.timehorizon=freepar(OCMATLSC.truncationtimecoord);
out.modelparameter=OCMATLSC.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='limitextremal';
out.solverinfo.truncationtimecoord=OCMATLSC.truncationtimecoord;
out.solverinfo.inftimetransformation=OCMATLSC.inftimetransformation;
out.solverinfo.pathtype=OCMATLSC.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATLSC.switchtimecoord;
if OCMATLSC.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATLSC.objectivevaluecoord;
end
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.solverinfo.yp=out.yp;
        out.solverinfo.ypmid=out.ypmid;
        out=rmfield(out,{'yp','ypmid'});
    case 'bvp5c'
    case 'bvp4c'
        out.solverinfo.yp=out.yp;
        out=rmfield(out,'yp');
end
if OCMATLSC.stateconstraint
    out.solverinfo.jumpcostatecoord=OCMATLSC.jumpcostatecoord;
end

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATLSC
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
%w=OCMATLSC.LAE_phi;
switch id
    case 1 % LP
        switch OCMATCONT.bvpmethod
            case 'sbvpoc'
            case 'bvp6c'
            case 'bvp4c'
                %                     [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                %                     dgdy=calchessianphi_bvp4c(tmesh,y,freepar,modelpar,odefun,bcfun,odejacfun,bcjacfun,@odehessian,@bchessian,v(1:OCMATCONT.HE.numdvariables-2),w(1:OCMATCONT.HE.numdvariables-2));
                %                     out(1)=dgdy(OCMATCONT.HE.numdvariables-2)'*w(OCMATCONT.HE.numdvariables-2);;
        end
        [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
        J=calc_RHSJac(t,y,z,freepar,modelpar,@ode,@bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,[]);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATLSC.PD_phi=p'/norm(p);
        OCMATLSC.PD_psi=Q(:,end);
        s.data.phi=OCMATLSC.PD_phi(:);
        s.data.laecoefficient=[];%nf_LAE(tmesh,coeff); % quadratic coefficient of center manifold
        s.data.sol=formatsolution(tmesh,coeff,tangent);
        s.msg =sprintf('Limit asymptotic extremal');
end
failed=0;

%------------------------------------------------------------
function [S,L]=singmat

% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

S=[0];

% S=[  0 8 8
%     8 0 8
%     1 8 0 ];

%L=[ 'BP'; 'H '; 'LP' ];
L=[ 'CP' ];


%elseif strcmp(arg, 'locate')
%--------------------------------------------------------
function [tmesh,coeff,tangent]=locate(id,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2)
switch id
    case 0
        [tmesh,coeff,tangent]=locateHP(tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2);
    otherwise
        error('No locator defined for singularity %d', id);
end
%---------------------------------------------------------
function done
WorkspaceDone;
% ---------------------------------------------------------
function WorkspaceDone(tmesh,coeff,tangent,varargin)

function workspaceadapt(tmesh,coeff,tangent)
% ------------------------------------------------------
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfunc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
J(:,OCMATCONT.HE.numdvariables+(-1:0))=[]; % remove derivative with respect to continuation parameter
J(OCMATCONT.HE.numdvariables-1,:)=[]; % remove derivative with respect to continuation parameter

[Q R E]=qr(full(J));
R1=R(1:end-1,1:end-1);
b=R(1:end-1,end);
w=E*[(R1\-b);1];
w=w/norm(w);
v=Q(:,end);

OCMATLSC.LAE_phi=w';
OCMATLSC.LAE_psi=v';
OCMATLSC.LAE_switch=1;



function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATLSC OCBVP

modelpar=OCMATLSC.parametervalue;
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        numarc=OCMATCONT.HE.numarc;
        domainddata=OCMATCONT.DOMAINDDATA;
        for arc=1:numarc
            arcindex=OCMATCONT.HE.arcindex(arc);
            idx=(OCMATCONT.HE.TIMEDDATA.leftarcindex(arc)-arc+1):(OCMATCONT.HE.TIMEDDATA.rightarcindex(arc)-arc);
            y(1:domainddata(arcindex).numode,1,idx)=coeff(OCMATCONT.HE.DDATA(arc).meshvalcoord); % values at the mesh points
            z(domainddata(arcindex).eqcoord,domainddata(arcindex).numcolscoord,idx)=coeff(OCMATCONT.HE.DDATA(arc).collvalcoord); % values at the collocation points
        end
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        % handle case of pure state constraints
    otherwise
end


function varargout=probleminit(varargin)
tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;

% ---------------------------------------------------------

function WorkspaceInit(tmesh,coeff,tangent)
global OCBVP

OCBVP.odehess=@odehess;
OCBVP.bchess=@bchess;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATLSC
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATLSC=OCMATLSC;
try
    if contnum==1
        save([OCMATLSC.basicglobalvarfilename '4limitextremal2lc'],'MODELINFO')
    end
    save([OCMATLSC.basicresultfilename '4limitextremal2lc'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATLSC

pathname=OCMATLSC.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
global OCMATLSC
% calculate phi and psi for next point
if OCMATLSC.LAE_switch == 0
    OCMATLSC.LAE_phi = OCMATLSC.LAE_new_phi/norm(OCMATLSC.LAE_new_phi);
else
    OCMATLSC.LAE_psi = OCMATLSC.LAE_new_psi/norm(OCMATLSC.LAE_new_psi);
end
OCMATLSC.LAE_switch = 1-OCMATLSC.LAE_switch;
flag=0;

function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATLSC
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATLSC.parametervalue);
        nN=n*N;
        % under the giving assumptions the solution is continuous even at
        % the border of the arcs. Therefore, the correspondent warning of
        % deval is turned 'off'
        warning('off','MATLAB:deval:NonuniqueSolution');
        ynew=devalbvpoc(sol,tmeshnew);
        warning('on','MATLAB:deval:NonuniqueSolution');
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,n,N,nN);

        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATLSC.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end
