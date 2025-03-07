function out=extremal2ep()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
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
% J=Jnum;
%J=full(J);
%if max(abs(Jnum(:)-J(:)))>1e-2
%    max(abs(Jnum(:)-J(:)))
%end

function opt=options(opt)

function varargout=probleminit(varargin)

tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(varargin);

% all done succesfully
varargout{1} = 0;

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
elseif OCMATAE.freeendtime
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATAE.objectivevaluecalc
    dxdt(OCMATAE.objectivevaluecoord,:)=dtds*OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
end
if OCMATAE.exogenousfunction
    dxdt(OCMATAE.exogenousdynamicscoordinate,:)=dtds*OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
elseif OCMATAE.freeendtime
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATAE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATAE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
if OCMATAE.exogenousfunction
    J=[[J zeros(size(J,1),1)];dtds*OCMATAE.exogenousjacobian(t,depvar,modelpar,arcarg)];
end

Jpar=zeros(OCMATAE.ODEcoord(end),OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATAE.autonomous
        Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATAE.objectivevaluecalc
        dxdt(OCMATAE.objectivevaluecoord,:)=OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATAE.objectivevaluecoord,:)=OCMATAE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATAE.exogenousfunction
        dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATAE.ODEcoord,OCMATAE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-(arc-1))*Jt;
        if arc>1
            Jpar(OCMATAE.ODEcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        Jpar(OCMATAE.ODEcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        if OCMATAE.movinghorizon
            dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
            if OCMATAE.objectivevaluecalc
                dxdt(OCMATAE.objectivevaluecoord,:)=OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
            end
            if OCMATAE.exogenousfunction
                dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
            end
            Jpar(:,OCMATAE.movinghorizoncoord)=dxdt;
        end
    end
else
    if OCMATAE.movinghorizon
        dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
        if OCMATAE.exogenousfunction
            dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
        end
        Jpar(:,OCMATAE.movinghorizoncoord)=dxdt;
    end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [Hy2,Hypar,Hpar2]=odehess(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
elseif OCMATAE.freeendtime
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
Hypar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
Hpar2=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter,OCMATCONT.HE.numparameter);

Hy2=dtds*OCMATAE.canonicalsystemhessian(t,depvar,modelpar,arcarg);
if OCMATCONT.HE.numarc>1
    J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
    if ~OCMATAE.autonomous
        Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if arc<OCMATCONT.HE.numarc
        Hypar(:,:,OCMATAE.switchtimecoord(arc))=J+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Hypar(:,:,OCMATAE.switchtimecoord(arc-1))=-(J+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        Hypar(:,:,OCMATAE.switchtimecoord(arc-1))=-(J+diffarctime(arc)*(s-arc)*Jt);
    end
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
switchtimes=freepar(OCMATAE.switchtimecoord);
if OCMATAE.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATAE.jumpcostateindex)=freepar(OCMATAE.jumpcostatecoord);
end
resconnec=[];
userbc=[];

if isfield(OCMATAE,'freevector') && ~isempty(OCMATAE.freevector)
    freevectorparameter=freepar(OCMATAE.freevectorindex);
    initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
    for ii=1:length(freevectorparameter)
        initialstate=initialstate+freevectorparameter(ii)*OCMATAE.freevector(:,ii);
    end
    
else
    initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
end
resinit=OCMATAE.bcinitial(depvara,OCMATAE.targetcoordinate,initialstate,modelpar,OCMATCONT.HE.arcarg(1));
try
    resasym=OCMATAE.bcasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint);
catch
    resasym=OCMATAE.bcasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint,modelpar,OCMATCONT.HE.arcarg(end),depvara);
end
if ~isempty(OCMATAE.fixedcoord)
    resasym=[resasym;depvarb(OCMATAE.fixedcoord,end)-OCMATAE.saddlepoint(OCMATAE.fixedcoord)];
end
if OCMATAE.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end
if OCMATAE.exogenousfunction
    resinit=[resinit; ...
        depvara(OCMATAE.exogenousdynamicscoordinate,1)-OCMATAE.exogenousinitialstates];
end
if OCMATAE.userbc
    if ~OCMATAE.movinghorizon
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
    elseif OCMATAE.freeendtime
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
    end
    userbc=OCMATAE.userbcfunc(arctime,depvara,depvarb,modelpar,OCMATCONT.HE.arcarg);
end

if OCMATAE.movinghorizon
    hatx=OCMATAE.saddlepoint;
    yend=depvarb(1:length(hatx),end);
    if ~isempty(OCMATAE.excludecoordinate4ep)
        yend(OCMATAE.excludecoordinate4ep,:)=[];
        hatx(OCMATAE.excludecoordinate4ep,:)=[];
    end
    resasym=[resasym; ...
        sqrt(sum((yend-hatx).^2))-OCMATAE.distance];
end
if OCMATAE.stateconstraint
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
else
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
end
res=[resinit;resconnec;resasym;userbc];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar);
    Ja(:,OCBVP.nBCs-OCBVP.nparmc+1:OCBVP.nBCs)=[];
    Jb(:,OCBVP.nBCs-OCBVP.nparmc+1:OCBVP.nBCs)=[];
    OCBVP.multiarccalc=0;
    return
end
switchtimes=freepar(OCMATAE.switchtimecoord);
if OCMATAE.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATAE.jumpcostateindex)=freepar(OCMATAE.jumpcostatecoord);
end
domainddata=OCMATCONT.DOMAINDDATA;

Jpar=zeros(OCMATCONT.HE.totalnumboundarycondition,OCMATCONT.HE.numparameter);
Ja=zeros(OCMATCONT.HE.totalnumboundarycondition-OCMATCONT.HE.numparameter+OCMATCONT.codimension);
Jb=Ja;
colcounter=0;
rowcounter=0;
arcindex=OCMATCONT.HE.arcindex(1);
colidx_start=colcounter+1;
colcounter=colcounter+domainddata(arcindex).numeq;
[Japart Jbpart Jpartmp]=OCMATAE.bcjacobianinitial(depvara,freepar,OCMATCONT.HE.arcarg(1),OCMATAE.targetcoordinate,OCMATAE.continuationvector);
if OCMATAE.objectivevaluecalc
    Japart(OCMATAE.objectivevaluecoord,OCMATAE.objectivevaluecoord)=1;
end

Jpar(1:OCMATCONT.HE.numinitialcondition,:)=Jpartmp;
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numinitialcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;
for arc=1:OCMATCONT.HE.numarc-1
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    Jaresetpart=OCMATAE.jacobianreset(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,-1);
    Jbresetpart=OCMATAE.jacobianreset(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,0);
    Jaguardpart=OCMATAE.jacobianguard(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,-1);
    Jbguardpart=OCMATAE.jacobianguard(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,0);
    
    % Ja
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+domainddata(arcindex).numode;
    Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbresetpart;
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+1;
    Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbguardpart;
    
    colidx_start=colcounter+1;
    colcounter=colcounter+domainddata(arcindex).numeq;
    % Jb
    rowcounter=rowcounter-domainddata(arcindex).numode-1; % 
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+domainddata(arcindex).numode;
    Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Jaresetpart;
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+1;
    Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Jaguardpart;
end

[Japart Jbpart  Jpartmp]=OCMATAE.bcjacobianasymptotic(depvarb,OCMATCONT.HE.arcarg(OCMATCONT.HE.numarc),OCMATAE.asymptoticmatrix',OCMATAE.saddlepoint);
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numendcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;
%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if ~OCMATAE.movinghorizon
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
    elseif OCMATAE.freeendtime
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1)); 
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATAE.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        %eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATAE.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    if any(strfind(OCMATAE.pathtype,'u'))
        violationmat=-diffarctime(arc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    else
        violationmat=diffarctime(arc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    end
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
if ~OCMATAE.movinghorizon
    yend=y(:,end);
    hatx=OCMATAE.saddlepoint;
    yend=yend(1:length(hatx));
    if ~isempty(OCMATAE.excludecoordinate4ep)
        yend(OCMATAE.excludecoordinate4ep,:)=[];
        hatx(OCMATAE.excludecoordinate4ep,:)=[];
    end
    %violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend(OCMATCONT.DOMAINDDATA(1).odecoord)-hatx)<0;
    violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx)<0;
    if violationmat
        counter=counter+1;
        cols=size(y,2);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows='maxdistance';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=norm(yend-hatx);
        infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx);
        b=min([b infoS(counter).minval]);
    end
end

%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
%sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
%figure(1)
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
elseif OCMATAE.freeendtime
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
for arc=1:numel(OCMATCONT.HE.arcarg)
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(arc)-diffarctime(arc)*(arc-1));
end

h=OCMATAE.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
%drawnow
%figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATAE
idx=[];
if isempty(coeff)
    return
end
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
if OCMATAE.stopcriterion
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    stopval=tangent(end);
    fprintf(1,' Stop value: %g\n',stopval);
    %fprintf(1,' Stop value: %g\n',OCMATAE.stopcriterionfunc(t,y,modelpar,[]));
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATAE
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
global OCMATCONT OCMATAE
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
            if isempty(tangent)
                return
            end
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
global OCMATCONT OCMATAE

failed=[];
for ii=id
    switch ii
        case 1
            if ~OCMATAE.stopcriterion
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
            else
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                out=tangent(end);
                %out=OCMATAE.stopcriterionfunc(t,y,modelpar,[]);
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATAE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
if ~OCMATAE.movinghorizon
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
elseif OCMATAE.freeendtime
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freeendtimecoord)];
    out.solverinfo.freeendtimecoord=OCMATAE.freeendtimecoord;
else
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
    out.solverinfo.distance=OCMATAE.distance;
    out.solverinfo.movinghorizoncoord=OCMATAE.movinghorizoncoord;
    %arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.switchtimecoord(end))+OCMATAE.movinghorizon];
end
%out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
if ~OCMATAE.movinghorizon
    out.timehorizon=OCMATAE.truncationtime;
else
    out.timehorizon=freepar(OCMATAE.movinghorizoncoord);
end
out.modelparameter=OCMATAE.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremal2ep';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.continuationvector=OCMATAE.continuationvector;
out.solverinfo.startvalue=OCMATAE.startvalue;
out.solverinfo.inftimetransformation=OCMATAE.inftimetransformation;
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATAE.switchtimecoord;

out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
if OCMATAE.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATAE.objectivevaluecoord;
end
if isfield(OCMATAE,'freevector') && ~isempty(OCMATAE.freevector)
    out.solverinfo.freevector=OCMATAE.freevector;
    out.solverinfo.freevectorindex=OCMATAE.freevectorindex;
else
    out.solverinfo.freevector=[];
    out.solverinfo.freevectorindex=[];
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
if OCMATAE.stateconstraint
    out.solverinfo.jumpcostatecoord=OCMATAE.jumpcostatecoord;
end

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        % reduce Jacobian to size without the continuation parameter
        J(:,OCMATCONT.HE.contparametercoord)=[];
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        q=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        q=q/norm(q);
        p=Q(:,end);
        p=p/norm(p);
        s.data.phi=q(:);
        s.data.psi=p(:);
        s.data.DFDX=J;
        try
            s.data.laecoefficient=calchnf_LP(tmesh,coeff,tangent,@operatoreq,q,p,J);
        catch
            s.data.laecoefficient=[];
        end
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
global OCMATCONT OCMATAE OCBVP
try
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    OCMATAE.probleminit(t,y,modelpar);
end
% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATAE OCBVP

modelpar=OCMATAE.parametervalue;
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


%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATAE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATAE=OCMATAE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATAE.basicglobalvarfilename '4extremal2ep'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremal2ep'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATAE

pathname=OCMATAE.datapath();

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATAE

discretizationdata=OCMATAE.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

