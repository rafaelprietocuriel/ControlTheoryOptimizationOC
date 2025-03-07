function out=modelperiodicsol()
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
% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end


function opt=options(opt)


function [opt out]=probleminit(tmesh,coeff,tangent)
global OCMATCONT OCMATPS
opt=setocoptions('OCCONTARG','WorkSpace',1);
opt=setocoptions(opt,'OCCONTARG','Locators',[1 0]);
out=0;


function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});


%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATPS OCMATCONT OCBVP
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
modelpar(OCMATPS.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
if OCMATPS.varyperiod
    period=OCMATPS.periodfunc(modelpar);
else
    period=OCMATPS.period;
end
arctime=[OCMATPS.initialtime freepar(OCMATPS.switchtimecoord).' OCMATPS.initialtime+period];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATPS.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATPS.linearizationcalc
    dtds*OCMATPS.canonicalsystemjacobian(t,depvar,modelpar,arcarg)*depvar(OCMATPS.linearizationindex);
    dxdt=[dxdt;ans(:)];
end
if OCMATPS.objectivevaluecalc
    dxdt(OCMATPS.objectivevaluecoord,:)=dtds*OCMATPS.objectivefunction(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATPS OCMATCONT OCBVP
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
modelpar(OCMATPS.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
if OCMATPS.varyperiod
    period=OCMATPS.periodfunc(modelpar);
else
    period=OCMATPS.period;
end
arctime=[OCMATPS.initialtime freepar(OCMATPS.switchtimecoord).' OCMATPS.initialtime+period];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATPS.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATPS.objectivevaluecalc
    J=[J; ...
        dtds*OCMATPS.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATPS.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATPS.autonomous
        Jt=OCMATPS.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATPS.objectivevaluecalc
        dxdt(OCMATPS.objectivevaluecoord,:)=OCMATPS.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATPS.objectivevaluecoord,:)=OCMATPS.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATPS.switchtimecoord(arc))=dxdt+dtds*(s-arc+1)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATPS.switchtimecoord(arc-1))=-(dxdt+dtds*(s-arc)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATPS.switchtimecoord(arc-1))=-(dxdt+dtds*(s-arc)*Jt);
    end
end
Jmodelpar=dtds*OCMATPS.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATPS.varyparameterindex);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATPS OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    res=bc(depvara,depvarb,freepar,modelpar);
    OCBVP.multiarccalc=0;
    return
end
modelpar(OCMATPS.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);

resconnec=[];
resinit=[];
resper=OCMATPS.bcperiodic(depvara,depvarb);
if OCMATPS.linearizationcalc
    resinit=depvara(OCMATPS.linearizationindex(:))-OCMATPS.linearizationinit;
end
if OCMATPS.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATPS.reset(depvara,depvarb,modelpar,freepar,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATPS.guard(depvara,depvarb,modelpar,freepar(OCMATPS.switchtimecoord),OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
end
res=[resconnec;resper;resinit];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,contval,arc)
global OCMATPS OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;
%----------------------------------------------------------------

function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATPS OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if OCMATPS.varyperiod
        period=OCMATPS.periodfunc(modelpar);
    else
        period=OCMATPS.period;
    end
    arctime=[OCMATPS.initialtime freepar(OCMATPS.switchtimecoord).' OCMATPS.initialtime+period];
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATPS.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATPS.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
global OCMATPS OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATPS.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
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
global OCMATCONT OCMATPS
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
global OCMATCONT OCMATPS
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
global OCMATCONT OCMATPS
failed=[];
out=[];
if isempty(coeff)
    failed=1;
    return
end

for ii=id
    switch ii
        case 1
            out=OCMATPS.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATPS OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out.octrajectory=transform2nativematlab(tmesh,coeff,OCMATPS.parametervalue);
if OCMATPS.linearizationcalc
    out.octrajectory.linearization=out.octrajectory.y(OCMATPS.linearizationindex(:),:);
    out.octrajectory.y(OCMATPS.linearizationindex(:),:)=[];
end
if OCBVP.multiarccalc
    out.octrajectory.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.octrajectory.arcposition=[1;length(out.octrajectory.x)];
end
if OCMATPS.varyperiod
    period=OCMATPS.periodfunc(modelpar);
else
    period=OCMATPS.period;
end
out.octrajectory.arcinterval=[OCMATPS.initialtime freepar(OCMATPS.switchtimecoord).' OCMATPS.initialtime+period];
out.octrajectory.arcarg=OCMATCONT.HE.arcarg;
out.octrajectory.x0=OCMATPS.initialtime;
out.octrajectory.timehorizon=OCMATPS.period;
out.period=OCMATPS.period;
out.octrajectory.modelparameter=OCMATPS.parametervalue;
out.octrajectory.modelparameter(OCMATPS.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);
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
out.octrajectory.linearization=calc_monodromy(t,y,z,freepar,modelpar,OCMATCONT.ode,OCMATCONT.odejac);
% add solver method specific information to make it consistent with MATLAB
% syntax
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATPS
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
% ------------------------------------------------------

function WorkspaceDone


%-----------------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATPS

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATPS.parametervalue;
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

modelpar(OCMATPS.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATPS

flag = 1;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATPS OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATPS=OCMATPS;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATPS.basicglobalvarfilename '4modelperiodicsol'],'MODELINFO')
    end
    save([OCMATPS.basicresultfilename '4modelperiodicsol'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATPS

discretizationdata=OCMATPS.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATPS

pathname=OCMATPS.datapath();


%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATPS
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATPS.parametervalue);
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
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATPS.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.solverinfo.yp;
        sol.idata.ymid=sol.solverinfo.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end

%-----------------------------------------------------------------
function M=monodromy(tmesh,coeff)
global OCMATCONT OCBVP

rows = OCBVP.rows;
cols = OCBVP.cols;
DPHI=OCMATCONT.frechetder(tmesh,coeff,[],OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac);
Jpi=DPHI(rows,cols);
cols=cols+OCBVP.numode;
Jpip1=DPHI(rows,cols);
M=-inv(Jpip1)*Jpi;
for ii=2:length(tmesh)-1
    rows = rows+OCBVP.numode;   % next equation
    Jpi=DPHI(rows,cols);
    cols=cols+OCBVP.numode;
    Jpip1=DPHI(rows,cols);
    M=-M*inv(Jpip1)*Jpi;
end
M=full(M);
