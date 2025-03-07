function out=indifferencedistribution()

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
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,icfunc},numel(coeff),numJacOpt);
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
global OCMATINDIF OCMATCONT OCBVP
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
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if OCMATINDIF.inftimetransformation(solutionindex) && arc==OCMATINDIF.cumsumnumarc(solutionindex)
    %t=arctime(arc)+log(arc-s)/OCMATINDIF.inftimetransformation;
    %dtds=-1./(OCMATINDIF.inftimetransformation*(arc-s));
else
    if solutionindex>1
        %timeshift=OCMATINDIF.truncationtime(solutionindex-1);
        transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
    else
        %timeshift=0;
        transformedtimeshift=0;
    end
    %t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc+1)-diffarctime(relarc)*relarc)+timeshift;
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    
    dtds=diffarctime(relarc);
end
dxdt=dtds*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
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
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
if OCMATINDIF.inftimetransformation(solutionindex) && arc==OCMATINDIF.cumsumnumarc(solutionindex)
    %t=arctime(arc)+log(arc-s)/OCMATINDIF.inftimetransformation;
    %dtds=-1./(OCMATINDIF.inftimetransformation*(arc-s));
else
    if solutionindex>1
        timeshift=OCMATINDIF.truncationtime(solutionindex-1);
        transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
    else
        timeshift=0;
        transformedtimeshift=0;
    end
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    dtds=diffarctime(relarc(ones(1,numel(s))));
end
J=OCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATINDIF.objectivevaluecalc
    J=[J; ...
        dtds*OCMATINDIF.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATINDIF.numarc(solutionindex)>1
    dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATINDIF.autonomous
        Jt=OCMATINDIF.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc))=dxdt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-dxdt;
        end
    else
        if OCMATINDIF.inftimetransformation
            %Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=-1./((exp(OCMATINDIF.inftimetransformation*arctime(arc))*(arc-s)))*dxdt;
        else
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        end
    end
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    res=bc(depvara,depvarb,freepar,modelpar);
    OCBVP.multiarccalc=0;
    return
end
resinit=[];
resasym=[];
resconnec=[];
residpt=[];

restarget=OCMATINDIF.bcinitial(depvara,OCMATINDIF.statecoordinate,OCMATINDIF.startvalue+freepar(end)*OCMATINDIF.continuationvector+freepar(OCMATINDIF.basisvectorcoord)*OCMATINDIF.basisvector);
for order=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.objectivevaluecalc
        resinit=[resinit;depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{order}(1))];
    end
    switchtimes=freepar(OCMATINDIF.switchtimecoord{order});
    if OCMATINDIF.stateconstraint
        jump=zeros(1,OCMATINDIF.numarc(order));
        jump(OCMATINDIF.jumpcostateindex{order})=freepar(OCMATINDIF.jumpcostatecoord{order});
    end
    arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{order});
    if order<OCMATINDIF.indifferenceorder
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(order+1))-depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(order))];
        residpt=[residpt; ...
            OCMATINDIF.bcindifference(depvara,modelpar,OCMATCONT.HE.arcarg([OCMATINDIF.arccoord{order:order+1}]),OCMATINDIF.initcoord([order order+1]))];
    end
    resasym=[resasym;OCMATINDIF.bcasymptotic(depvarb(:,OCMATINDIF.cumsumnumarc(order)),OCMATINDIF.asymptoticmatrix{order},OCMATINDIF.saddlepoint{order})];
    if OCMATINDIF.stateconstraint
        for arc=1:OCMATINDIF.numarc(order)-1
            resconnec=[resconnec; ...
                OCMATINDIF.reset(depvara(:,OCMATINDIF.arccoord{order}),depvarb(:,OCMATINDIF.arccoord{order}),modelpar,switchtimes,jump,arcarg,OCMATINDIF.edge{order},arc); ...
                OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{order}),depvarb(:,OCMATINDIF.arccoord{order}),modelpar,switchtimes,jump,arcarg,OCMATINDIF.edge{order},arc)];
        end
    else
        for arc=1:OCMATINDIF.numarc(order)-1
            resconnec=[resconnec; ...
                OCMATINDIF.reset(depvara(:,OCMATINDIF.arccoord{order}),depvarb(:,OCMATINDIF.arccoord{order}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{order},arc); ...
                OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{order}),depvarb(:,OCMATINDIF.arccoord{order}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{order},arc)];
            if OCMATINDIF.objectivevaluecalc
                resconnec=[resconnec;depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{order}(arc+1))-depvarb(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{order}(arc))];
            end
        end
    end
end
res=[restarget;resinit;resasym;resconnec;residpt];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
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

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATINDIF OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATINDIF.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));    
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATINDIF.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATINDIF.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    violationmat=diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    if any(violationmat)
        counter=counter+1;
        cols=find(violationmat);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows='switchtime';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        %infoS(counter).constraintvalue=diffarctime(arc);
        infoS(counter).minval=min(diffarctime(:));
        b=min([b infoS(counter).minval]);
    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
if length(sol.x)~=size(sol.y,2)
    sol
end
% clear possible persistent variable
h=OCMATINDIF.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
%global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
fprintf(1,' Continuation parameter: %g\n',coeff(end));

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATINDIF
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
global OCMATCONT OCMATINDIF
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
global OCMATCONT

failed=[];
for ii=id
    switch ii
        case 1
            out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,OCMATINDIF.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{1}).' OCMATINDIF.truncationtime(1), ...
    OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{2}).' OCMATINDIF.truncationtime(2)];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATINDIF.initialtime;
if OCMATINDIF.inftimetransformation(1)
    out.timehorizon=inf;
else
    out.timehorizon=OCMATINDIF.truncationtime;
end
out.modelparameter=OCMATINDIF.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencedistribution';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATINDIF.inftimetransformation;
out.solverinfo.pathtype=OCMATINDIF.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.basisvectorcoord=OCMATINDIF.basisvectorcoord;
out.solverinfo.basisvector=OCMATINDIF.basisvector;
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
if OCMATINDIF.stateconstraint
    out.solverinfo.jumpcostatecoord=OCMATINDIF.jumpcostatecoord;
end


%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATINDIF
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATINDIF.PD_phi=p'/norm(p);
        OCMATINDIF.PD_psi=Q(:,end);
        s.data.phi=OCMATINDIF.PD_phi(:);
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
global OCMATCONT OCMATINDIF

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATINDIF.parametervalue;
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

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATINDIF OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATINDIF=OCMATINDIF;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATINDIF.basicglobalvarfilename '4indifferencedistribution'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencedistribution'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATINDIF

discretizationdata=OCMATINDIF.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATINDIF

pathname=OCMATINDIF.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATINDIF
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATINDIF.parametervalue);
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
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATINDIF.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end