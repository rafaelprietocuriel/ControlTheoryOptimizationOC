function out=limitcycle()
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
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
res=calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
% J=full(J);

% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,icfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=Jnum;
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end

function opt=options()
global OCMATLC
opt=setocoptions('OCCONTARG','WorkSpace',1,'Adapt',1);
if isfield(OCMATLC,'targetparametervalue') && isempty(OCMATLC.targetparametervalue)
    opt=setocoptions(opt,'OCCONTARG','HitTargetValue',0);
end
opt=opt.MATCONT;  

function varargout=probleminit(varargin)
tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;



function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});


%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
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
modelpar(OCMATLC.varyparameterindex)=freepar(OCMATLC.parametercoord);
arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATLC.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
function out=icfun(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    fun=zeros(size(depvar));
    OCBVP.multiarccalc=1;
    for ii=1:OCMATCONT.HE.numarc
        depvarint=depvar(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:);
        fun(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:)=icfun(s,depvarint,ii,freepar,modelpar);
    end
    OCBVP.multiarccalc=0;
    return
end
%arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
%diffarctime=diff(arctime);
%dtds=diffarctime(arc);

if ~isempty(OCMATLC.oldlc.solver)
    [dum fun]=devalbvpoc(OCMATLC.oldlc,s);
else
    fun=interp1(OCMATLC.oldlc.x,OCMATLC.oldlc.yp.',s).';
end
out=dot(fun,depvar);
% 
%-------------------------------------------------------------------------
function [J Jpar]=icfunjac(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    fun=zeros(size(depvar));
    OCBVP.multiarccalc=1;
    for ii=1:OCMATCONT.HE.numarc
        depvarint=depvar(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:);
        fun(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:)=icfun(s,depvarint,ii,freepar,modelpar);
    end
    OCBVP.multiarccalc=0;
    return
end
%arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
%diffarctime=diff(arctime);
%dtds=diffarctime(arc);

if ~isempty(OCMATLC.oldlc.solver)
    [dum fun]=devalbvpoc(OCMATLC.oldlc,s);
else
    fun=interp1(OCMATLC.oldlc.x,OCMATLC.oldlc.yp.',s).';
end
J=fun.';

Jpar=zeros(OCBVP.numic,OCBVP.npar);
%Jpar(1,end-1)=sum(fun.*depvar);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATLC.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATLC.autonomous
        Jt=OCMATLC.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
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


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
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
arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).'  freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATLC.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATLC.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATLC.autonomous
        Jt=OCMATLC.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
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
Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATLC.varyparameterindex);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    res=bc(depvara,depvarb,freepar,modelpar);
    OCBVP.multiarccalc=0;
    return
end
modelpar(OCMATLC.varyparameterindex)=freepar(OCMATLC.parametercoord);
switchtimes=freepar(OCMATLC.switchtimecoord);
if OCMATLC.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATLC.jumpcostateindex)=freepar(OCMATLC.jumpcostatecoord);
end

resconnec=[];
resinit=[];
resper=OCMATLC.bcperiodic(depvara,depvarb);
if OCMATLC.stateconstraint
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATLC.reset(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATLC.guard(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
else
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATLC.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATLC.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
end
res=[resconnec;resper;resinit];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
[Ja Jb Jpar]=OCMATLC.bcjacobianperiodic(depvara,depvarb);


%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;
%----------------------------------------------------------------

function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATLC OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATLC.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)
%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));

%---------------------------------------------------------------
function varargout=defaultprocessor(tmesh,coeff,tangent)
global OCMATCONT OCMATLC
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
modelpar(OCMATLC.varyparameterindex)=freepar(OCMATLC.parametercoord);
arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);

%dataadaptation(tmesh);
for arc=1:OCMATCONT.HE.numarc
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    arcarg=OCMATCONT.HE.arcarg(arc);
    depvarint=y(OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc),:);
    OCMATLC.oldlc.yp(OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc),:)=OCMATLC.canonicalsystem(t,depvarint,modelpar,arcarg);
    %OCMATLC.oldlc.yp(OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc),:)=freepar(OCMATLC.periodcoord)*OCMATLC.canonicalsystem(t,depvarint,modelpar,arcarg);
end
OCMATLC.oldlc.parameters=freepar(OCMATLC.periodcoord);
OCMATLC.oldlc.x=s;
OCMATLC.oldlc.y=y;


varargout{2}=nan;
% all done succesfully
varargout{1}=0;
%-------------------------------------------------------------
function [out, failed]=testfunc(id,tmesh,coeff,tangent)
global OCMATCONT OCMATLC
%
%if any(ismember(id,[1 2]))
%    J=cjac(cds.curve_func,cds.curve_jacobian,x,[]);
%end

out(1)=0;
failed=[];

for ii=id
    lastwarn('');

    switch ii
        %case 1 % BP
        % Jalcobian extended with bordering vectors v and w
        %B=[J; v'];
        %out(1)=det(B);


        case 1 % LP
            out(1)=tangent(end);

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
            if ~isempty(OCMATLC.targetparametervalue)
                out=OCMATLC.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATLC OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out.octrajectory=transform2nativematlab(tmesh,coeff,OCMATLC.parametervalue);
if OCMATLC.linearizationcalc
    out.octrajectory.linearization=out.octrajectory.y(OCMATLC.linearizationindex(:),:);
    out.octrajectory.y(OCMATLC.linearizationindex(:),:)=[];
end
if OCBVP.multiarccalc
    out.octrajectory.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.octrajectory.arcposition=[1;length(out.octrajectory.x)];
end
out.octrajectory.arcinterval=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).'  freepar(OCMATLC.periodcoord)];
out.octrajectory.arcarg=OCMATCONT.HE.arcarg;
out.octrajectory.x0=OCMATLC.initialtime;
out.octrajectory.timehorizon=freepar(OCMATLC.periodcoord);
out.period=freepar(OCMATLC.periodcoord);
out.octrajectory.modelparameter=OCMATLC.parametervalue;
out.octrajectory.modelparameter(OCMATLC.varyparameterindex)=freepar(OCMATLC.parametercoord);
out.octrajectory.modelname=OCMATCONT.modelname;

out.octrajectory.solver=OCMATCONT.bvpmethod;
out.octrajectory.solverinfo.coeff=coeff;
out.octrajectory.solverinfo.tmesh=tmesh;
out.octrajectory.solverinfo.tangent=tangent;
out.octrajectory.solverinfo.parameters=freepar;
out.octrajectory.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.octrajectory.solverinfo.multiarccalc=OCBVP.multiarccalc;
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.octrajectory.solverinfo.yp=out.octrajectory.yp;
        out.octrajectory.solverinfo.ypmid=out.octrajectory.ypmid;
        out.octrajectory=rmfield(out.octrajectory,{'yp','ypmid'});
    case 'bvp5c'
    case 'bvp4c'
        out.octrajectory.solverinfo.yp=out.octrajectory.yp;
        out.octrajectory=rmfield(out.octrajectory,'yp');
end
if OCMATLC.stateconstraint
    out.octrajectory.solverinfo.jumpcostatecoord=OCMATLC.jumpcostatecoord;
end
out.octrajectory.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.octrajectory.solverinfo.switchtimecoord=OCMATLC.switchtimecoord;
out.octrajectory.solverinfo.periodcoord=OCMATLC.periodcoord;
out.octrajectory.solverinfo.oldlc=OCMATLC.oldlc;
out.octrajectory.solverinfo.oldlcInit=OCMATLC.oldlcInit;
out.octrajectory.linearization=calc_monodromy(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.odejac);
% add solver method specific information to make it consistent with MATLAB
% syntax
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATLC
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
global OCMATCONT OCMATLC
OCMATLC.oldlc=OCMATLC.oldlcInit;
% ------------------------------------------------------

function WorkspaceDone


%-----------------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATLC

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATLC.parametervalue;
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
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
    otherwise
end

modelpar(OCMATLC.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATLC
flag = 1;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATLC OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATLC=OCMATLC;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATLC.basicglobalvarfilename '4limitcycle'],'MODELINFO')
    end
    save([OCMATLC.basicresultfilename '4limitcycle'],'sout','bvpout')
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
