function out=heteroclinicep2fep()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;
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
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,icfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end

function opt=options(opt)

function varargout=probleminit(varargin)
tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP

solutionindex=OCMATHET.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if ~OCMATHET.freeendtime(solutionindex)
    arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
else
    arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.truncationtimecoord{solutionindex}).'];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATHET.numarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
solutionindex=OCMATHET.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if ~OCMATHET.freeendtime(solutionindex)
    arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
else
    arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.truncationtimecoord{solutionindex}).'];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
if solutionindex>1
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);
else
    timeshift=0;
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATHET.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATHET.numarc(solutionindex)>1
    dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATHET.autonomous
        Jt=OCMATHET.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if relarc<OCMATHET.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoord{solutionindex}(relarc))=dxdt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoord{solutionindex}(relarc-1))=-dxdt;
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
    end
end
if OCMATHET.freeendtime(solutionindex)
    dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
    Jpar(:,OCMATHET.truncationtimecoord{solutionindex})=dxdt;
end
Jmodelpar=dtds*OCMATHET.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATHET.parametervaluecoord)=Jmodelpar(:,OCMATHET.parameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP

resinit=[];
resasym=[];
resconnec=[];
resequilibrium=[];
algebraicequation=[];

modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

for ii=1:OCMATHET.hetorder
    switchtimes=freepar(OCMATHET.switchtimecoord{ii});
    arcarg=OCMATCONT.HE.arcarg(OCMATHET.arccoord{ii});
    if ii<OCMATHET.hetorder
        resinit=[resinit; ...
            depvara(OCMATHET.statecostatecoordinate,OCMATHET.initcoord(ii+1))-depvara(OCMATHET.statecostatecoordinate,OCMATHET.initcoord(ii))-freepar(end)*OCMATHET.initialstatedifference];
        if ~isempty(OCMATHET.fixcoordinate)
            resinit=[resinit; ...
                OCMATHET.fixvalue-depvara(OCMATHET.fixcoordinate,OCMATHET.initcoord(ii))];
        end
    end
    if ~isempty(OCMATHET.algebraicequation)
        algebraicequation=[algebraicequation; ...
            OCMATHET.algebraicequation(depvara(:,OCMATHET.initcoord(ii)),modelpar,arcarg(1))];
    end
    hatx=freepar(OCMATCONT.HE.equilibriumcoord{ii});
    if ~OCMATHET.ftequilibriumindex(ii)
        resequilibrium=[resequilibrium;OCMATHET.equilibrium(hatx,modelpar,arcarg(end))];
    else
        resequilibrium=[resequilibrium;OCMATHET.explicitequilibriumvalue(hatx,modelpar,arcarg(end))];
    end
    if ~isempty(OCMATHET.algebraicequation)
        resequilibrium=[resequilibrium; ...
            OCMATHET.algebraicequation(hatx,modelpar,arcarg(end))];
    end
    resasym=[resasym;OCMATHET.bcasymptotic(depvarb(:,OCMATHET.cumsumnumarc(ii)),OCMATHET.asymptoticmatrix{ii},hatx)];
    for arc=1:OCMATHET.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATHET.reset(depvara(:,OCMATHET.arccoord{ii}),depvarb(:,OCMATHET.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc); ...
            OCMATHET.guard(depvara(:,OCMATHET.arccoord{ii}),depvarb(:,OCMATHET.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc)];
    end
end
res=[resinit;resasym;resconnec;resequilibrium;algebraicequation];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

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
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp6c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp5c(t,y,tangent,res,canRemovePoints,freepar,modelpar,ode);
        coeffnew=[ynew(:)];
    case 'bvp4c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp4c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
end

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATHET OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATHET.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATHET.numarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    if ~OCMATHET.freeendtime(solutionindex)
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    else
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.truncationtimecoord{solutionindex}).'];
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATHET.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATHET.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
%     if OCMATHET.stableflag{arc}
%         violationmat=diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
%     else
%         violationmat=-diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
%     end
%     if any(violationmat)
%         counter=counter+1;
%         cols=find(violationmat);
%         infoS(counter).arcarg=arcarg;
%         infoS(counter).arcnum=arc;
%         infoS(counter).rows='switchtime';
%         infoS(counter).cols=cols;
%         infoS(counter).violationmat=violationmat;
%         %infoS(counter).constraintvalue=diffarctime(arc);
%         infoS(counter).minval=min(diffarctime(:));
%         b=min([b infoS(counter).minval]);
%     end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
if length(sol.x)~=size(sol.y,2)
    sol
end
% clear possible persistent variable
h=OCMATHET.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
%global OCMATHET OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
fprintf(1,' Continuation parameter: %g\n',coeff(end));

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET

failed=[];
for ii=id
    switch ii
        case 1
            out=OCMATHET.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
            %out=tangent(end);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATHET OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,OCMATHET.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[];
for ii=1:OCMATHET.hetorder
    if ~OCMATHET.freeendtime(ii)
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{ii}).' OCMATHET.truncationtime(ii)];
    else
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{ii}).' freepar(OCMATHET.truncationtimecoord{ii}).'];
    end
    out.arcinterval=[out.arcinterval arctime];
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATHET.initialtime;
out.timehorizon=OCMATHET.truncationtime;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='heteroclinicep2fep';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATHET.inftimetransformation;
out.solverinfo.pathtype=OCMATHET.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATHET.switchtimecoord;
out.solverinfo.truncationtimecoord=OCMATHET.truncationtimecoord;
out.solverinfo.equilibriumcoord=OCMATCONT.HE.equilibriumcoord;
out.solverinfo.parameterindex=OCMATHET.parameterindex;
out.solverinfo.parametervaluecoord=OCMATHET.parametervaluecoord;
out.solverinfo.ftequilibriumindex=OCMATHET.ftequilibriumindex;
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.solverinfo.yp=out.yp;
        out.solverinfo.ypmid=out.ypmid;
        out=rmfield(out,{'yp','ypmid'});
    case 'bvp4c'
        out.solverinfo.yp=out.yp;
        out=rmfield(out,'yp');
end


%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATHET
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
s=[];
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


% ------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATHET

%numarc=OCMATCONT.HE.numarc;
%domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATHET.parametervalue;
if length(OCMATHET.parameterindex)==2
    modelpar(OCMATHET.parameterindex)=coeff(OCMATCONT.HE.contparametercoord+(-1:0));
else
    modelpar(OCMATHET.parameterindex)=coeff(OCMATCONT.HE.contparametercoord-1);
end
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATHET OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATHET=OCMATHET;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATHET.basicglobalvarfilename '4heteroclinicep2fep'],'MODELINFO')
    end
    save([OCMATHET.basicresultfilename '4heteroclinicep2fep'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATHET

discretizationdata=OCMATHET.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATHET

pathname=OCMATHET.datapath();


%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATHET

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
%
for ii=1:OCMATHET.hetorder
    arcarg=OCMATCONT.HE.arcarg(OCMATHET.arccoord{ii});
    hatx=freepar(OCMATCONT.HE.equilibriumcoord{ii});
    if ~OCMATHET.ftequilibriumindex(ii)
        J{ii}=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
        OCMATHET.asymptoticmatrix{ii}=asymptoticbc(J{ii},OCMATHET.pathtype{ii},'c');
    end
end
flag = 1;

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATHET
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATHET.parametervalue);
        nN=n*N;
        % under the giving assumptions the solution is continuous even at
        % the border of the arcs. Therefore, the correspondent warning of
        % deval is turned 'off'
        warning('off','MATLAB:deval:NonuniqueSolution');
        ynew=devalbvpoc(sol,tmeshnew);
        warning('on','MATLAB:deval:NonuniqueSolution');
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,n,N,nN);

        coeffnew=[ynew(:);freepar(:)];
end


%-----------------------------------------------------------------
function out=ricatti(Y,J,Q0,dim)
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(Q0,J,dim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
