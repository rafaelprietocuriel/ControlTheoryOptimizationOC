function out=limitlimitcycle()
% continuation file for a periodic solution of a non-autonomous system

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{4}{3}=@icfun;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;
out{5}{3}=@icfunjac;

out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
out{11}=@operatorpfrechet;
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{17}=@adapt;
out{18}=@dataadaptation;
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
out{25}=@deval;
out{26}=@saveintermediate;
out{27}=@datapath;
out{28}=@domaindiscretization;
out{30}=@printcontinuation;

function res=operatoreq(tmesh,coeff,tangent,odefun,bcfun,icfunc)
global OCMATCONT OCMATLC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
res=calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
% append g (minimally augmented system, e.g. Kuzentsov 1998, p. 502f, 2004, p. 544f (3rd edition))
M=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfunc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
M=reduceJac(M,OCMATLC.LPC_phi,OCMATLC.LPC_psi);
%M=[J  psi]
%  [phi' 0]
% the function g is defined as
%       M*[w;g]=[0;1], w\in\R^N, (*)
% subsequently we set s=[w;g] and b=[0;1].
% If s is a solution of eq. (*) then J*w+g*phi=0. Thus, for g=0, J*w=0 and w is
% a right eigenvector of J for the eigenvalue zero.
% If s is a solution of the transposed eq. (*) then v*J'+g*psi=0, with
% v=w'. Thus, for g=0 v is a left eigenvector of J for the eigenvalue zero.

% update phi and psi alternating to keep M non-singular
% psi not in Im M, phi not in Im M'
% these properties are guaranteed by the following procedure
b = []; b(OCMATCONT.HE.numdvariables-1)=1;
if OCMATLC.LPC_switch==1
    % solves M'*s=b => s in kern M'= (Im M)^+
    s = M'\b';
    OCMATLC.LPC_new_psi = s(1:OCMATCONT.HE.numdvariables-2)';
else
    % if J has a single zero eigenvector then J* has a also a single zero
    % eigenvector, therefore M'*s=b' does the same job as M*s=b'
    % solves M*s=b => s in Im M= (Kern M)^+
    s = M\b';
    OCMATLC.LPC_new_phi= s(1:OCMATCONT.HE.numdvariables-2)';
end
res(OCMATCONT.HE.numdvariables-1,1)=s(OCMATCONT.HE.numdvariables-1);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)
global OCMATCONT OCMATLC

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);

M=reduceJac(J,OCMATLC.LPC_phi,OCMATLC.LPC_psi);

% calculate g'
b=[]; 
b(OCMATCONT.HE.numdvariables-1)=1;
w=M\b';
v=M'\b';

switch OCMATCONT.bvpmethod
    case 'bvp4c'
        dgdy=calcgprime_bvp4c(tmesh,y,freepar,modelpar,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.odehess,OCMATCONT.bchess,v,w);
end

J(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables)=dgdy;

function opt=options(opt)
global OCMATLC
opt=setocoptions(opt,'OCCONTARG','WorkSpace',1,'Adapt',1);
if isempty(OCMATLC.targetparametervalue)
    opt=setocoptions(opt,'OCCONTARG','HitTargetValue',0);
end

function varargout=probleminit(varargin)
tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;



function [F,J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});


%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
modelpar(OCMATLC.freeparameterindex)=freepar(OCMATLC.freeparametercoord);
arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATLC.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATLC.objectivevaluecalc
    dxdt(OCMATLC.objectivevaluecoord,:)=dtds*OCMATLC.objectivefunction(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; derivatives of the form dtds*f_t have to be included

function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
modelpar(OCMATLC.freeparameterindex)=freepar(OCMATLC.freeparametercoord);

arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).'  freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATLC.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
dxdt=OCMATLC.canonicalsystem(t,depvar,modelpar,arcarg);
if ~OCMATLC.autonomous
    Jt=OCMATLC.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
else
    Jt=0;
end
if OCMATLC.objectivevaluecalc
    dxdt(OCMATLC.objectivevaluecoord,:)=OCMATLC.objectivefunction(t,depvar,modelpar,arcarg);
    if ~OCMATLC.autonomous
        Jt(OCMATLC.objectivevaluecoord,:)=OCMATLC.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
end
if OCMATCONT.HE.numarc==1
    Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLC.periodcoord)=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
else
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLC.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLC.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLC.periodcoord)=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLC.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
    end
end
Jmodelpar=dtds*OCMATLC.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATLC.freeparametercoord)=Jmodelpar(:,OCMATLC.freeparameterindex);
%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
modelpar(OCMATLC.freeparameterindex)=freepar(OCMATLC.freeparametercoord);
switchtimes=freepar(OCMATLC.switchtimecoord);

resconnec=[];
%
resper=OCMATLC.bcperiodic(depvara,depvarb);
if isempty(OCMATLC.fixcoordinate)
    resper=[resper; ...
    sum(OCMATLC.velocityvector.*(depvara(OCMATLC.velocitycoord,1)-OCMATLC.initialpoint))];
else
    resper=[resper; ...
    depvara(OCMATLC.fixcoordinate,1)-OCMATLC.fixvalue];
end
if OCMATLC.objectivevaluecalc
    resper=[resper;depvara(OCMATLC.objectivevaluecoord,1)];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATLC.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATLC.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
end
res=[resconnec;resper];

%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
global OCMATCONT OCMATLC OCBVP
failed=0;
N=length(tmeshold);
neqnN=OCBVP.numode*length(tmeshold);
%[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
mbcidx = find(diff(tmeshold) == 0);  % locate internal interfaces
ismbvp = ~isempty(mbcidx);
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(tmeshold)];
nregions = length(mbcidx) + 1;
mbcidxint = find(diff(tmesh) == 0);  % locate internal interfaces
Lidxint = [1, mbcidxint+1];
Ridxint = [mbcidxint, length(tmesh)];
phifreepar=OCMATLC.LPC_phi(neqnN+1:end);
phi=reshape(OCMATLC.LPC_phi(1:neqnN),OCBVP.numode,N);
phiint=[];
psifreepar=OCMATLC.LPC_psi(neqnN+1:end);
psi=reshape(OCMATLC.LPC_psi(1:neqnN),OCBVP.numode,N);
psiint=[];
if ismbvp
    for region=1:nregions
        xidx=Lidx(region):Ridx(region);
        xidxint=Lidxint(region):Ridxint(region);
        phiint=[phiint interp1(tmeshold(xidx),phi(:,xidx).',tmesh(xidxint)).'];
        psiint=[psiint interp1(tmeshold(xidx),psi(:,xidx).',tmesh(xidxint)).'];
    end
else
    phiint=interp1(tmeshold,phi.',tmesh).';
    psiint=interp1(tmeshold,psi.',tmesh).';
end
phi=[phiint(:);phifreepar(:)];
psi=[psiint(:);psifreepar(:)];

%phi=interp1(tmeshold,OCMATLC.LPC_phi(OCMATCONT.HE.DDATA.meshvalcoordold)',tmesh).';
%psi=interp1(tmeshold,OCMATLC.LPC_psi(OCMATCONT.HE.DDATA.meshvalcoordold)',tmesh).';
phi=phi/norm(phi);
psi=psi/norm(psi);
% calc_RHS(t,y,z,freepar,modelpar,@ode,@bc,[]);
% %append g (minimally augmented system, e.g. Kuzentsov 1998, p. 502f)
% M=calc_RHSJac(t,y,z,freepar,modelpar,@ode,@bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
% M=reduceJac(M,phi(:),psi(:));
% 
% b = []; b(OCMATCONT.HE.numdvariables-1,1)=1;
% phi = M\b;
% phi=phi(1:OCMATCONT.HE.numdvariables-2)';
% phi=phi/norm(phi);
% psi = M'\b;
% psi=psi(1:OCMATCONT.HE.numdvariables-2)';
% psi=psi/norm(psi);

OCMATLC.LPC_phi_old=OCMATLC.LPC_phi;
OCMATLC.LPC_psi_old=OCMATLC.LPC_psi;
OCMATLC.LPC_phi=phi(:);
OCMATLC.LPC_psi=psi(:);
OCMATLC.LPC_new_phi=phi(:);
OCMATLC.LPC_new_psi=psi(:);
%----------------------------------------------------------------

function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATLC OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATLC.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATLC.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
    end
    violationmat=constr<-OCMATCONT.OPTIONS.admissibletol;

    negarctime=diff(arctime)<-OCMATCONT.OPTIONS.admissibletol;
    if any(violationmat(:)) || any(negarctime)
        counter=counter+1;
        [rows cols]=find(violationmat);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows=rows;
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=constr;
        infoS(counter).arctime=arctime;
        infoS(counter).minval=min(constr(:));
        b=min([b infoS(counter).minval diff(arctime)]);
    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATLC OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
%sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
%figure(1)
s=s(leftarcindex(1):rightarcindex(end));
t=s;
y=y(:,leftarcindex(1):rightarcindex(end));
arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
for arc=1:OCMATCONT.HE.numarc
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*(s(leftarcindex(arc):rightarcindex(arc)))+(arctime(arc)-diffarctime(arc)*(arc-1));
end

% clear possible persistent variable
%figure(1)
h=OCMATLC.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
% drawnow
% figure(gcf)
%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATLC
idx=[];
if isempty(coeff)
    return
end
if isempty(OCMATLC.targetparametervalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(end));
else
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
    fprintf(1,' Difference to targetvalue: %g\n',OCMATLC.targetparametervalue-modelpar(OCMATLC.targetparameterindex));
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATLC
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
global OCMATCONT OCMATLC

out(1)=0;
failed=[];

for ii=id
    lastwarn('');

    switch ii
        case 1 % LP
            switch OCMATCONT.bvpmethod
                case 'sbvpoc'
                case 'bvp6c'
                case 'bvp4c'
                    out(1)=calchnf_LP(tmesh,coeff,tangent,@operatoreq,OCMATLC.LPC_phi,OCMATLC.LPC_psi);
                    OCMATLC.testval=out;
            end

        otherwise
            error('No such testfunction');
    end
    if ~isempty(lastwarn)
        failed=[failed ii];
    end

end
%----------------------------------------------------------------
function [out, failed]=targetvaluefunc(id,tmesh,coeff,tangent)
global OCMATCONT OCMATLC
failed=[];
out=[];
if isempty(coeff)
    failed=1;
    return
end

for ii=id
    switch ii
        case 1
            if isempty(OCMATLC.targetparametervalue)
                out=coeff(end);
            else
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                out=OCMATLC.targetparametervalue-modelpar(OCMATLC.targetparameterindex);
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATLC OCBVP
%dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out.octrajectory=transform2nativematlab(tmesh,coeff,OCMATLC.parametervalue);
if OCBVP.multiarccalc
    out.octrajectory.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.octrajectory.arcposition=OCMATCONT.HE.arcrowindex;
end
out.octrajectory.arcinterval=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).'  freepar(OCMATLC.periodcoord)];
out.octrajectory.arcarg=OCMATCONT.HE.arcarg;
out.octrajectory.x0=OCMATLC.initialtime;
out.octrajectory.timehorizon=freepar(OCMATLC.periodcoord);
out.period=freepar(OCMATLC.periodcoord);
out.octrajectory.modelparameter=modelpar;
out.octrajectory.modelname=OCMATCONT.modelname;

out.octrajectory.solver=OCMATCONT.bvpmethod;
out.octrajectory.solverinfo.coeff=coeff;
out.octrajectory.solverinfo.tmesh=tmesh;
out.octrajectory.solverinfo.tangent=tangent;
out.octrajectory.solverinfo.parameters=freepar;
out.octrajectory.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
switch OCMATCONT.bvpmethod
    case 'bvp4c'
        out.octrajectory.solverinfo.yp=out.octrajectory.yp;
        out.octrajectory=rmfield(out.octrajectory,'yp');
end
out.octrajectory.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.octrajectory.solverinfo.switchtimecoord=OCMATLC.switchtimecoord;
out.octrajectory.solverinfo.periodcoord=OCMATLC.periodcoord;
if OCMATLC.monodromy
    out.octrajectory.linearization=calc_monodromy(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.odejac);
    out.octrajectory.solverinfo.Hi=OCBVP.Hi;
    out.octrajectory.solverinfo.Gi=OCBVP.Gi;
else
    out.octrajectory.linearization=[];
end
% add solver method specific information to make it consistent with MATLAB
% syntax
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATLC
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
failed=0;
%s.data=[];
s.msg =sprintf('Limit asymptotic extremal');
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
L=[ 'LP' ];


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

function WorkspaceInit(tmesh,coeff,tangent)
global OCMATCONT OCMATLC OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

calc_RHS(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
J=calc_RHSJac(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
J(:,OCMATCONT.HE.numdvariables+(-1:0))=[]; % remove derivative with respect to continuation parameter
J(OCMATCONT.HE.numdvariables-1,:)=[]; % remove derivative with respect to continuation parameter

[Q R E]=qr(full(J));
R1=R(1:end-1,1:end-1);
b=R(1:end-1,end);
w=E*[(R1\-b);1];
w=w/norm(w);
v=Q(:,end);

OCMATLC.LPC_phi=v';
OCMATLC.LPC_psi=w';
OCMATLC.LPC_switch=1;
% ------------------------------------------------------

function WorkspaceDone


%-----------------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATLC

modelpar=OCMATLC.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end

modelpar(OCMATLC.freeparameterindex)=freepar(OCMATLC.freeparametercoord);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATLC 
%[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

%OCMATLC.initialpoint=coeff(OCMATLC.velocitycoord);
%OCMATLC.velocityvector=ode(0,y(:,1),1,freepar,modelpar);
%OCMATLC.velocityvector=OCMATLC.velocityvector/norm(OCMATLC.velocityvector);

% calculate phi and psi for next point
if OCMATLC.LPC_switch == 0
    OCMATLC.LPC_phi = OCMATLC.LPC_new_phi/norm(OCMATLC.LPC_new_phi);
else
    OCMATLC.LPC_psi = OCMATLC.LPC_new_psi/norm(OCMATLC.LPC_new_psi);
end
OCMATLC.LPC_switch = 1-OCMATLC.LPC_switch;
flag=0;


%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATLC OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATLC=OCMATLC;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATLC.basicglobalvarfilename '4limitlimitcycle'],'MODELINFO')
    end
    save([OCMATLC.basicresultfilename '4limitlimitcycle'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATLC

discretizationdata=OCMATLC.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATLC

pathname=OCMATLC.datapath();


%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATLC
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATLC.parametervalue);
        nN=n*N;
        % under the giving assumptions the solution is continuous even at
        % the border of the arcs. Therefore, the correspondent warning of
        % deval is turned 'off'
        warning('off','MATLAB:deval:NonuniqueSolution');
%         sol.yp=sol.solverinfo.yp;
%         try
%             sol.ypmid=sol.solverinfo.ypmid;
%         end
        ynew=devalbvpoc(sol,tmeshnew);
        warning('on','MATLAB:deval:NonuniqueSolution');
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,n,N,nN);

        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATLC.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.solverinfo.yp;
        sol.idata.ymid=sol.solverinfo.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end


function J=reduceJac(J,phi,psi)
global OCMATCONT

J(:,OCMATCONT.HE.numdvariables)=[]; % remove derivative with respect to continuation parameter
J(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=psi(:);
J(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=phi(:).';
J(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;


