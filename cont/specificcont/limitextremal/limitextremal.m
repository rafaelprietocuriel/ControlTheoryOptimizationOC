function out=limitextremal()

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
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{17}=@adapt;
out{18}=@dataadaptation;
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
% append g (minimally augmented system, e.g. Kuzentsov 1998, p. 502f, 2004, p. 544f (3rd edition))
M=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfunc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
M=reduceJac(M,OCMATLSC.LAE_phi,OCMATLSC.LAE_psi);
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
if OCMATLSC.LAE_switch==1
    % solves M'*s=b => s in kern M'= (Im M)^+
    s = M'\b';
    OCMATLSC.LAE_new_psi = s(1:OCMATCONT.HE.numdvariables-2)';
else
    % if J has a single zero eigenvector then J* has a also a single zero
    % eigenvector, therefore M'*s=b' does the same job as M*s=b'
    % solves M*s=b => s in Im M= (Kern M)^+
    s = M\b';
    OCMATLSC.LAE_new_phi= s(1:OCMATCONT.HE.numdvariables-2)';
end
res(OCMATCONT.HE.numdvariables-1,1)=s(OCMATCONT.HE.numdvariables-1);
res(OCMATCONT.HE.numdvariables,1)=0;

function res=operatoreqTest(tmesh,coeff,tangent,odefun,bcfun,icfunc)
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
%res=calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
% append g (minimally augmented system, e.g. Kuzentsov 1998, p. 502f, 2004, p. 544f (3rd edition))
M=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfunc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
M=reduceJac(M,OCMATLSC.LAE_phi,OCMATLSC.LAE_psi);
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
s = M\b';
res=s(OCMATCONT.HE.numdvariables-1);
%res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
M=reduceJac(J,OCMATLSC.LAE_phi,OCMATLSC.LAE_psi);

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

% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% dy=1e-6;
% for ii=1:length(coeff)
%     coeff0=coeff;
%     coeff0(ii)=coeff0(ii)+dy;
%     resh=operatoreqTest(tmesh,coeff0,tangent,odefun,bcfun,icfun);
%     coeff0=coeff;
%     coeff0(ii)=coeff0(ii)-dy;
%     resl=operatoreqTest(tmesh,coeff0,tangent,odefun,bcfun,icfun);
%     dgdynum(ii)=(resh-resl)/2/dy;
% end
% dgdynum
%Jnum=numjaccsd(@operatoreqTest,{tmesh,coeff,tangent,odefun,bcfun,icfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end



function J=reduceJac(J,phi,psi)
global OCMATCONT

J(:,OCMATCONT.HE.numdvariables)=[]; % remove derivative with respect to continuation parameter
J(1:OCMATCONT.HE.numdvariables-2,OCMATCONT.HE.numdvariables-1)=psi(:);
J(OCMATCONT.HE.numdvariables-1,1:OCMATCONT.HE.numdvariables-2)=phi(:).';
J(OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables-1)=0;


function [F,J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
if ~OCMATLSC.freeendtime
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
else
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.freeendtimecoord)];
end
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
% has to be adapted;  derivatives of the form dtds*f_t have to be included

function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP

if ~OCMATLSC.freeendtime
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
else
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.freeendtimecoord)];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATLSC.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATLSC.objectivevaluecalc
    J=[J; ...
        dtds*OCMATLSC.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATLSC.canonicalsystem(t,depvar,modelpar,arcarg);
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
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLSC.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLSC.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATLSC.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
    end
end
if OCMATLSC.freeendtime
    dxdt=OCMATLSC.canonicalsystem(t,depvar,modelpar,arcarg);
    Jpar(:,OCMATLSC.freeendtimecoord)=dxdt;
end


%-------------------------------------------------------------------------
function [Hy2,Hypar,Hpar2]=odehess(s,depvar,arc,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
if ~OCMATLSC.freeendtime
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
else
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.freeendtimecoord)];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
Hypar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
Hpar2=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter,OCMATCONT.HE.numparameter);

Hy2=dtds*OCMATLSC.canonicalsystemhessian(t,depvar,modelpar,arcarg);
if OCMATCONT.HE.numarc>1
    J=OCMATLSC.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
    if ~OCMATLSC.autonomous
        Jt=OCMATLSC.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if arc<OCMATCONT.HE.numarc
        Hypar(:,:,OCMATLSC.switchtimecoord(arc))=J+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Hypar(:,:,OCMATLSC.switchtimecoord(arc-1))=-(J+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        Hypar(:,:,OCMATLSC.switchtimecoord(arc-1))=-(J+diffarctime(arc)*(s-arc)*Jt);
    end
end
if OCMATLSC.freeendtime
    J=OCMATLSC.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
    Hypar(:,:,OCMATLSC.freeendtimecoord)=J;
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
switchtimes=freepar(OCMATLSC.switchtimecoord);
resconnec=[];

initialstate=OCMATLSC.startvalue;
for ii=1:length(OCMATLSC.freevectorcoord)
    initialstate=initialstate+freepar(OCMATLSC.freevectorcoord(ii))*OCMATLSC.freevector(:,ii);
end
resinit=depvara(OCMATLSC.statecoord,1)-initialstate;
resasym=OCMATLSC.bcasymptotic(depvarb,OCMATLSC.asymptoticmatrix,OCMATLSC.saddlepoint);
if OCMATLSC.freeendtime>0
    resasym=[resasym; ...
        sqrt(sum((OCMATLSC.saddlepoint-depvarb(:,end)).^2))-OCMATLSC.distance];
end
if OCMATLSC.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATLSC.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATLSC.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
end
res=[resinit;resconnec;resasym];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjacobian(depvara,depvarb,freepar,modelpar)
global OCMATLSC OCMATCONT
Ja=[];
Jb=[];
Jpar=[];


%-------------------------------------------------------------------------
function [Ha Hb Hpar]=bchess(depvara,depvarb,freepar,modelpar)
global OCMATLSC OCMATCONT OCBVP
Ha=zeros(OCBVP.nBCs,2*OCBVP.numode+OCBVP.nparmcod,OCBVP.numode);
Hb=Ha;
Hpar=zeros(OCBVP.nBCs,2*OCBVP.numode+OCBVP.nparmcod,OCBVP.npar);

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
global OCMATCONT OCMATLSC OCBVP
failed=0;
N=length(tmeshold);
neqnN=OCBVP.numode*length(tmeshold);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
mbcidx = find(diff(tmeshold) == 0);  % locate internal interfaces
ismbvp = ~isempty(mbcidx);
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(tmeshold)];
nregions = length(mbcidx) + 1;
mbcidxint = find(diff(tmesh) == 0);  % locate internal interfaces
Lidxint = [1, mbcidxint+1];
Ridxint = [mbcidxint, length(tmesh)];
phifreepar=OCMATLSC.LAE_phi(neqnN+1:end);
phi=reshape(OCMATLSC.LAE_phi(1:neqnN),OCBVP.numode,N);
phiint=[];
psifreepar=OCMATLSC.LAE_psi(neqnN+1:end);
psi=reshape(OCMATLSC.LAE_psi(1:neqnN),OCBVP.numode,N);
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

%phi=interp1(tmeshold,OCMATLSC.LAE_phi(OCMATCONT.HE.DDATA.meshvalcoordold)',tmesh).';
%psi=interp1(tmeshold,OCMATLSC.LAE_psi(OCMATCONT.HE.DDATA.meshvalcoordold)',tmesh).';
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

OCMATLSC.LAE_phi_old=OCMATLSC.LAE_phi;
OCMATLSC.LAE_psi_old=OCMATLSC.LAE_psi;
OCMATLSC.LAE_phi=phi(:);
OCMATLSC.LAE_psi=psi(:);
OCMATLSC.LAE_new_phi=phi(:);
OCMATLSC.LAE_new_psi=psi(:);


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
    if ~OCMATLSC.freeendtime
        arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
    else
        arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.freeendtimecoord)];
    end
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
global OCMATCONT OCMATLSC
idx=[];
if isempty(coeff)
    return
end

fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
try
    fprintf(1,'           NF parameter: %g\n',OCMATLSC.testval);
end
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

for ii=id
    lastwarn('');
    switch ii
        case 1 % ECP
            switch OCMATCONT.bvpmethod
                case 'sbvpoc'
                case 'bvp6c'
                case 'bvp4c'
                    out(1)=calchnf_LP(tmesh,coeff,tangent,@operatoreq,OCMATLSC.LAE_phi,OCMATLSC.LAE_psi);
                    OCMATLSC.testval=out;
            end
    end
    if 0%~isempty(lastwarn)
        failed=[failed ii];
    end

end
%----------------------------------------------------------------
function [out, failed]=targetvaluefunc(id,tmesh,coeff,tangent)
global OCMATCONT OCMATLSC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

failed=[];
for ii=id
    switch ii
        case 1
            if ~isempty(OCMATLSC.targetcoordinate)
                out=OCMATLSC.targetvalue-y(OCMATLSC.targetcoordinate,1);
            else
                out=[];
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATLSC OCBVP
%dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATLSC.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
if ~OCMATLSC.freeendtime
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' OCMATLSC.truncationtime];
else
    arctime=[OCMATLSC.initialtime freepar(OCMATLSC.switchtimecoord).' freepar(OCMATLSC.freeendtimecoord)];
end
out.arcinterval=arctime;
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATLSC.initialtime;
out.timehorizon=OCMATLSC.truncationtime;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='limitextremal';
%out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATLSC.inftimetransformation;
out.solverinfo.pathtype=OCMATLSC.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATLSC.switchtimecoord;
out.solverinfo.freevector=OCMATLSC.freevector;
out.solverinfo.freevectorcoord=OCMATLSC.freevectorcoord;
out.solverinfo.data.psi=OCMATLSC.LAE_psi;
out.solverinfo.data.phi=OCMATLSC.LAE_phi;

out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
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
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATLSC
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));

switch id
    case 1 % CP
        [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
        J=calc_RHSJac(t,y,z,freepar,modelpar,@ode,@bc,[],OCMATCONT.odejac,OCMATCONT.bcjac,[]);
        J1=J(1:end-1,1:end-2);
        [Q,R,E]=qr(full(J1)); % see Govaerts et al 2005 (p. 242)
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATLSC.PD_phi=p'/norm(p);
        OCMATLSC.PD_psi=Q(:,end);
        s.data.phi=OCMATLSC.PD_phi(:);
        s.data.psi=OCMATLSC.PD_psi(:);
        
        s.data.laecoefficient=calchnf_CP(tmesh,coeff,tangent,@operatoreq,OCMATLSC.PD_phi,OCMATLSC.PD_psi,J1);
        s.data.sol=formatsolution(tmesh,coeff,tangent);
        s.msg =sprintf('Cusp point');
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
global OCMATCONT OCMATLSC OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

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

OCMATLSC.LAE_phi=v';
OCMATLSC.LAE_psi=w';
OCMATLSC.LAE_switch=1;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATLSC OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATLSC=OCMATLSC;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATLSC.basicglobalvarfilename '4limitextremal'],'MODELINFO')
    end
    save([OCMATLSC.basicresultfilename '4limitextremal'],'sout','bvpout')
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
