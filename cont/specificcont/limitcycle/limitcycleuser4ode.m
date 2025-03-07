function out=limitcycleuser4ode()
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
out{21}=@plotcontinuation;
out{22}=@formatsolution;
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

% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,icfunc},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end

function opt=options(opt)
global OCMATLC
opt=setocoptions(opt,'OCCONTARG','WorkSpace',1,'Adapt',1,'CheckAdmissibility','off');
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



function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});


%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP

modelpar(OCMATLC.varyparameterindex)=freepar(OCMATLC.parametercoord);
arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).' freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATLC.dynamics(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP

arctime=[OCMATLC.initialtime freepar(OCMATLC.switchtimecoord).'  freepar(OCMATLC.periodcoord)];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATLC.dynamicsjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATLC.dynamics(t,depvar,modelpar,arcarg);
    if ~OCMATLC.autonomous
        Jt=OCMATLC.dynamicsderivativetime(t,depvar,modelpar,arcarg);
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
Jmodelpar=dtds*OCMATLC.dynamicsparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATLC.varyparameterindex);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
modelpar(OCMATLC.varyparameterindex)=freepar(OCMATLC.parametercoord);
switchtimes=freepar(OCMATLC.switchtimecoord);

resconnec=[];
resinit=[];
% resper=[OCMATLC.bcperiodic(depvara,depvarb); ...
%     depvara(OCMATLC.fixedcoord,1)-OCMATLC.fixedvalue];
resper=[OCMATLC.bcperiodic(depvara,depvarb); ...
    sum(OCMATLC.velocityvector.*(depvara(:,1)-OCMATLC.initialpoint))];
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATLC.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATLC.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
end
res=[resconnec;resper;resinit];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATLC OCMATCONT OCBVP
[Ja Jb Jpar]=OCMATLC.bcjacobianperiodic(depvara,depvarb);

%----------------------------------------------------------------
function [res,varargout]=residual(tmesh,coeff,rhs,odefun)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
OCMATCONT.interpolateintegralconstraintfunc=true;
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
OCMATCONT.interpolateintegralconstraintfunc=false;

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

switch OCMATCONT.bvpmethod
    case 'bvp5c'
        b=verifyresidual_bvp5c(maxres);
    otherwise
        b=maxres < OCMATCONT.OPTIONS.meshadaptreltol;
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

%----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=meshadaptation(tmesh,coeff,tangent,res,canRemovePoints,ode)
global OCMATCONT OCMATLC
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp4c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp4c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp5c(t,y,tangent,res,canRemovePoints,freepar,modelpar,ode);
        coeffnew=[ynew(:)];
    case 'bvp6c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp6c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
end
%----------------------------------------------------------------
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
            out=OCMATLC.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
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
out.octrajectory.arcposition=[1;length(out.octrajectory.x)];

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
out.octrajectory.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.octrajectory.solverinfo.switchtimecoord=OCMATLC.switchtimecoord;
out.octrajectory.solverinfo.periodcoord=OCMATLC.periodcoord;
%out.octrajectory.linearization=calc_monodromy(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.odejac);
out.octrajectory.linearization=[];
%out.octrajectory.solverinfo.Hi=OCBVP.Hi;
%out.octrajectory.solverinfo.Gi=OCBVP.Gi;
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
%global OCMATCONT OCMATLC
%OCMATLC.oldlc=OCMATLC.oldlcInit;
% ------------------------------------------------------

function WorkspaceDone


%-----------------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATLC

numarc=OCMATCONT.HE.numarc;
modelpar=OCMATLC.parametervalue;
switch OCMATCONT.bvpmethod
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
%dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
OCMATLC.initialpoint=coeff(OCMATLC.velocitycoord);
OCMATLC.velocityvector=ode(0,y(:,1),1,freepar,modelpar);
OCMATLC.velocityvector=OCMATLC.velocityvector/norm(OCMATLC.velocityvector);

%OCMATLC.oldlc=transform2nativematlab(tmesh,coeff,OCMATLC.parametervalue);
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
        save([OCMATLC.basicglobalvarfilename '4limitcycleuser'],'MODELINFO')
    end
    save([OCMATLC.basicresultfilename '4limitcycleuser'],'sout','bvpout')
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
