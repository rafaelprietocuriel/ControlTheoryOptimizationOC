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

function [F,J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP

switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.freeparameterindex)=freemodelparameter;
    otherwise
        if ~isempty(OCMATAE.freeparameter)
            modelpar(OCMATAE.freeparameterindex)=freepar(OCMATAE.freeparametercoordinate);
        end
end

transformedtimeshift=0;
if strcmp(OCMATCONT.continuationtype,'time')
    arctime=OCMATAE.switchtimes;
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    %arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
    if OCMATAE.freeendtime
        arctime=[0 arctime freepar(OCMATAE.endtimecoordinate)];
    else
        arctime=[0 arctime OCMATAE.endtime];
    end
else
    arctime=OCMATAE.switchtimes;
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    if OCMATAE.freeendtime
        arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
    else
        arctime=[0 arctime OCMATAE.endtime];
    end
end
arcarg=OCMATAE.arcarg(arc);
diffarctime=diff(arctime);
t=diffarctime(arc)*(s-transformedtimeshift)+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATAE.objectivevaluecalc
    %dxdt(OCMATAE.objectivevaluecoordinate(:,arc),:)=dtds*OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
    dxdt(OCMATAE.objectivevaluecoordinate,:)=dtds*OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
end
if OCMATAE.exogenousfunction
    dxdt(OCMATAE.exogenousdynamicscoordinate,:)=dtds*OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
transformedtimeshift=0;
if strcmp(OCMATCONT.continuationtype,'time')
    arctime=OCMATAE.switchtimes;
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    %arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
    if OCMATAE.freeendtime
        arctime=[0 arctime freepar(OCMATAE.endtimecoordinate)];
    else
        arctime=[0 arctime OCMATAE.endtime];
    end
else
    arctime=OCMATAE.switchtimes;
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    if OCMATAE.freeendtime
        arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
    else
        arctime=[0 arctime OCMATAE.endtime];
    end
end
arcarg=OCMATAE.arcarg(arc);
diffarctime=diff(arctime);
t=diffarctime(arc)*(s-transformedtimeshift)+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);

J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if length(OCBVP.numode)>1
    numode=OCBVP.numode(arc);
else
    numode=OCBVP.numode;
end
if OCMATAE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATAE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(numode,length(OCMATAE.objectivevaluecoordinate))];
end
if OCMATAE.exogenousfunction
    J=[[J zeros(size(J,1),OCMATAE.exogenousnumberofstates)];[dtds*OCMATAE.exogenousjacobian(t,depvar,modelpar,arcarg),zeros(OCMATAE.exogenousnumberofstates)]];
end
Jpar=zeros(numode,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATAE.autonomous
        Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATAE.objectivevaluecalc
        dxdt(OCMATAE.objectivevaluecoordinate,:)=OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATAE.objectivevaluecoordinate,:)=OCMATAE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATAE.exogenousfunction
        dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
    end
    if arc<OCMATAE.numarc
        if 1%any(arc==OCMATAE.freetimeindex)
            Jpar(1:numode,OCMATAE.switchtimecoordinate(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        end
        if arc>1
            if 1%any(arc-1==OCMATAE.freetimeindex)
                Jpar(1:numode,OCMATAE.switchtimecoordinate(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
            end
        end
    else
        if 1%any(arc==OCMATAE.freetimeindex)
            Jpar(1:numode,OCMATAE.switchtimecoordinate(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
        if OCMATAE.freeendtime
            Jpar(1:numode,OCMATAE.endtimecoordinate)=dxdt;
            if any(arc-1==OCMATAE.freetimeindex)
                %Jpar(1:numode,OCMATAE.switchtimecoordinate(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
            end
        end
    end
else
    if OCMATAE.freeendtime
        dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
        if OCMATAE.exogenousfunction
            dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
        end
        Jpar(:,OCMATAE.endtimecoordinate)=dxdt;
    end
end
switch OCMATCONT.continuationtype
    case 'parameter'
        if OCMATAE.implicit
            Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
            %Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg,OCMATAE.freeparameterindex);
        else
            Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
        end
        if isempty(OCMATAE.targetvalue)
            Jpar(:,OCMATAE.continuationcoordinate)=Jmodelpar(:,OCMATAE.freeparameterindex);
        else
            Jpar(OCMATAE.statecostatecoordinate,OCMATCONT.HE.numparameter)=OCMATAE.continuationvector.'*Jmodelpar(:,OCMATAE.continuationparameterindex);
        end
        if OCMATAE.exogenousfunction
            Jmodelpar=dtds*OCMATAE.exogenousparameterjacobian(t,depvar,modelpar,arcarg);
            Jpar(OCMATAE.exogenousdynamicscoordinate,OCMATCONT.HE.numparameter)=OCMATAE.continuationvector.'*Jmodelpar(:,OCMATAE.continuationparameterindex);
        end
    case 'time'
        if ~OCMATAE.autonomous
            Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
        else
            Jt=0;
        end
        if arc==OCMATAE.continuationtimeindex
            Jpar(1:numode,OCMATAE.continuationcoordinate)=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        elseif any(arc-1==OCMATAE.continuationtimeindex)
            Jpar(1:numode,OCMATAE.continuationcoordinate)=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
end
if ~isempty(OCMATAE.freeparameter) && ~strcmp(OCMATCONT.continuationtype,'parameter')
    if OCMATAE.implicit
        Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg,OCMATAE.freeparameterindex);
    else
        Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
    end
    Jpar(:,OCMATAE.continuationcoordinate)=Jmodelpar(:,OCMATAE.freeparameterindex);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [Hy2,Hypar,Hpar2]=odehess(s,depvar,arc,freepar,modelpar)
Hy2=[];
Hypar=[];
Hpar2=[];

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
resinit=[];
resconnec=[];
resricatti=[];
resequilibrium=[];
resuser=[];

switchtimes=OCMATAE.switchtimes;
switchtimes(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.freeparameterindex)=freemodelparameter;
        saddlepoint=freepar(OCMATAE.EP.equilibriumcoordinate);
        initialstate=OCMATAE.initialpoint(OCMATAE.statecoordinate);
        Jac=OCMATAE.canonicalsystemjacobian(0,saddlepoint,modelpar,OCMATAE.EP.arcarg);
        if OCMATAE.EP.implicitcontrolnum
            dudx=OCMATAE.dimplicitcontroldx(0,saddlepoint,modelpar,OCMATAE.EP.arcarg);
            Jac=Jac(OCMATAE.statecostatecoordinate,OCMATAE.statecostatecoordinate)+Jac(OCMATAE.statecostatecoordinate,OCMATAE.statecostatecoordinate(end)+(1:OCMATAE.EP.implicitcontrolnum))*dudx;
        end

        if ~OCMATAE.simple
            Y=freepar(OCMATAE.Ycoordinate);
            asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
            resricatti=ricatti(Y,Jac);
        else
            %asymptoticmatrix=asymptoticbc(Jac,OCMATAE.pathtype,'c');
            asymptoticmatrix=OCMATAE.asymptoticmatrix;
        end
        resequilibrium=OCMATAE.equilibrium(saddlepoint,modelpar,OCMATAE.EP.arcarg);
        initialstatecoordinate=OCMATAE.initialstatecoordinate;

    case 'initialstate'
        initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
        if isempty(OCMATAE.freeparameter)
            asymptoticmatrix=OCMATAE.asymptoticmatrix;
        end
        saddlepoint=OCMATAE.EP.saddlepoint;
        initialstatecoordinate=OCMATAE.depvarcoordinate;
    case 'time'
        switchtimes(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        initialstate=OCMATAE.initialpoint(OCMATAE.statecoordinate);
        asymptoticmatrix=OCMATAE.asymptoticmatrix;
        saddlepoint=OCMATAE.EP.saddlepoint;
        initialstatecoordinate=OCMATAE.initialstatecoordinate;
end
if ~isempty(OCMATAE.freeparameter) && ~OCMATAE.simple
    Y=freepar(OCMATAE.Ycoordinate);
    asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
    saddlepoint=freepar(OCMATAE.EP.equilibriumcoordinate);
    Jac=OCMATAE.canonicalsystemjacobian(0,saddlepoint,modelpar,OCMATAE.EP.arcarg);
    if OCMATAE.EP.implicitcontrolnum
        dudx=OCMATAE.dimplicitcontroldx(0,saddlepoint,modelpar,OCMATAE.EP.arcarg);
        Jac=Jac(OCMATAE.statecostatecoordinate,OCMATAE.statecostatecoordinate)+Jac(OCMATAE.statecostatecoordinate,OCMATAE.statecostatecoordinate(end)+(1:OCMATAE.EP.implicitcontrolnum))*dudx;
    end
    resricatti=ricatti(Y,Jac);
    resequilibrium=OCMATAE.equilibrium(saddlepoint,modelpar,OCMATAE.EP.arcarg);
end

if ~isempty(OCMATAE.targetcontrolcoordinate)
    resinit=[resinit; ...
        OCMATAE.bccontrol(depvara(:,OCMATAE.trajectoryindex(1)),OCMATAE.targetcontrolcoordinate, ...
        OCMATAE.initialcontrolvalue+freepar(end)*OCMATAE.controlvaluedirection,modelpar,OCMATAE.arcarg(1))];
end
if ~isempty(OCMATAE.fixinitstate)
    %resinit=depvara(OCMATAE.fixinitstate,OCMATAE.trajectoryindex(1))-OCMATAE.initialpoint(OCMATAE.fixinitstate);
end
if ~isempty(initialstatecoordinate)
    resinit=[resinit; ...
        OCMATAE.bcinitial(depvara(:,OCMATAE.trajectoryindex(1)),initialstatecoordinate,initialstate,modelpar,OCMATAE.arcarg(1))];
end
arcarg=OCMATAE.arcarg;

resasym=OCMATAE.bcasymptotic(depvarb(:,OCMATAE.trajectoryindex(end)),asymptoticmatrix,saddlepoint);
if OCMATAE.fixdistance
    resasym=[resasym; ...
        sqrt(sum((saddlepoint(OCMATAE.fixdistancecoordinate)-depvarb(OCMATAE.fixdistancecoordinate,end)).^2))-OCMATAE.distance];
%     resasym=[resasym; ...
%         norm(saddlepoint(OCMATAE.fixdistancecoordinate)-depvarb(OCMATAE.fixdistancecoordinate,OCMATAE.trajectoryindex(end)))-OCMATAE.distance];
end
if OCMATAE.objectivevaluecalc
    resinit=[resinit;depvara(OCMATAE.objectivevaluecoordinate,1)];
end
if OCMATAE.exogenousfunction
    resinit=[resinit; ...
        depvara(OCMATAE.exogenousdynamicscoordinate,1)-OCMATAE.exogenousinitialstates];
end
if OCMATAE.userbc
    resuser=OCMATAE.userbcfunc(depvara,depvarb,modelpar,OCMATCONT.HE.arcarg);
end

for arc=1:OCMATAE.numarc-1
    if any(arc==OCMATAE.freetimeindex)
        resconnec=[resconnec; ...
            OCMATAE.reset(depvara(:,OCMATAE.trajectoryindex),depvarb(:,OCMATAE.trajectoryindex),modelpar,switchtimes,arcarg,OCMATAE.edge,arc); ...
            OCMATAE.guard(depvara(:,OCMATAE.trajectoryindex),depvarb(:,OCMATAE.trajectoryindex),modelpar,switchtimes,arcarg,OCMATAE.edge,arc)];
    else
        resconnec=[resconnec; ...
            OCMATAE.reset(depvara(:,OCMATAE.trajectoryindex),depvarb(:,OCMATAE.trajectoryindex),modelpar,switchtimes,arcarg,OCMATAE.edge,arc)];
    end
    if OCMATAE.objectivevaluecalc
        %resconnec=[resconnec;depvara(OCMATAE.objectivevaluecoordinate(arc+1),arc+1)-depvarb(OCMATAE.objectivevaluecoordinate(arc),arc)];
        resconnec=[resconnec;depvara(OCMATAE.objectivevaluecoordinate,arc+1)-depvarb(OCMATAE.objectivevaluecoordinate,arc)];
    end
    if OCMATAE.exogenousfunction
        resconnec=[resconnec;depvara(OCMATAE.exogenousdynamicscoordinate,arc+1)-depvarb(OCMATAE.exogenousdynamicscoordinate,arc)];
    end
    
end

res=[resinit;resconnec;resasym;resricatti;resequilibrium;resuser];

%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;
%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;

for arc=1:OCMATCONT.HE.numarc
    if strcmp(OCMATCONT.continuationtype,'time')
        arctime=OCMATAE.switchtimes;
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
        %arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        if OCMATAE.freeendtime
            arctime=[0 arctime freepar(OCMATAE.endtimecoordinate)];
        else
            arctime=[0 arctime OCMATAE.endtime];
        end
    else
        arctime=OCMATAE.switchtimes;
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
        if OCMATAE.freeendtime
            arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
        else
            arctime=[0 arctime OCMATAE.endtime];
        end
    end
    diffarctime=diff(arctime);
    arcarg=OCMATAE.arcarg(arc);
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
checkcoordinate=setdiff(OCMATAE.statecostatecoordinate,OCMATAE.fixdistancecoordinate);
if ~isempty(checkcoordinate) & ~strcmp(OCMATCONT.continuationtype,'time')
    switch OCMATCONT.continuationtype
        case 'parameter'
            saddlepoint=freepar(OCMATAE.EP.equilibriumcoordinate);
        otherwise
            if ~isempty(OCMATAE.freeparameter)
                saddlepoint=freepar(OCMATAE.EP.equilibriumcoordinate);
            else
                saddlepoint=OCMATAE.EP.saddlepoint;
            end
    end
    if OCMATAE.implicit
        saddlepoint=saddlepoint(OCMATAE.statecostatecoordinate);
    end
    yend=y(1:length(saddlepoint),end);

    violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend(checkcoordinate)-saddlepoint(checkcoordinate))<0;
    if violationmat
        counter=counter+1;
        cols=size(y,2);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows='maxdistance';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=norm(yend-saddlepoint);
        infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(yend-saddlepoint);
        b=min([b infoS(counter).minval]);
    end
end

%----------------------------------------   ------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
%sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
%figure(1)
s=s(leftarcindex(1):rightarcindex(end));
t=s;
y=y(:,leftarcindex(1):rightarcindex(end));
if strcmp(OCMATCONT.continuationtype,'time')
    arctime=[OCMATAE.switchtimes OCMATAE.endtime];
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    %arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
    arctime=[0 arctime];
    if OCMATAE.freeendtime
        arctime(end)=freepar(OCMATAE.endtimecoordinate);
    end
else
    arctime=OCMATAE.switchtimes;
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    if OCMATAE.freeendtime
        arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
    else
        arctime=[0 arctime OCMATAE.endtime];
    end
end
diffarctime=diff(arctime);
for arc=1:OCMATCONT.HE.numarc
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*(s(leftarcindex(arc):rightarcindex(arc)))+(arctime(arc)-diffarctime(arc)*(arc-1));
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
switch OCMATCONT.continuationtype
    case {'initialstate','parameter'}
        if ~OCMATAE.hitfunction
            fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
        else
            [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
            depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
            depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);

            out=OCMATAE.targetfunction(depvara,depvarb,modelpar,OCMATAE.arcarg);
            fprintf(1,' Hit Value: %g\n',out);
        end
    case 'time'
        if ~isempty(OCMATAE.targetvalue)
            [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
            depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
            depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
            switch OCMATAE.targettype
                case 'T'
                    out=OCMATAE.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
                case 'd'
                    out=OCMATAE.targetvalue-norm(depvarb(OCMATAE.statecostatecoordinate,end)-OCMATAE.EP.saddlepoint(OCMATAE.statecostatecoordinate));
            end
            fprintf(1,' Hit Value: %g\n',out);

        end

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
function [out,failed]=targetvaluefunc(id,tmesh,coeff,tangent)
global OCMATCONT OCMATAE

failed=[];
for ii=id
    switch ii
        case 1
            switch OCMATCONT.continuationtype
                case {'initialstate','parameter'}
                    if ~OCMATAE.hitfunction
                        out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
                    else
                        [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                        depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                        depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                        out=OCMATAE.targetfunction(depvara,depvarb,modelpar,OCMATAE.arcarg);
                    end
                case 'time'
                    if ~isempty(OCMATAE.targetvalue)
                        [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                        depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                        depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                        switch OCMATAE.targettype
                            case 'T'
                                out=OCMATAE.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
                            case 'd'
                                [t,y,z,freepar]=drearr(tmesh,coeff);
                                out=OCMATAE.targetvalue-norm(depvarb(OCMATAE.statecostatecoordinate,end)-OCMATAE.EP.saddlepoint(OCMATAE.statecostatecoordinate));
                        end

                    else
                        out=[];
                    end
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
switch OCMATCONT.continuationtype
    case 'initialstate'
        arctime=OCMATAE.switchtimes;
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
        if OCMATAE.freeendtime
            arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
        else
            arctime=[0 arctime OCMATAE.endtime];
        end
        timehorizon=arctime(end);
    case 'parameter'
        arctime=OCMATAE.switchtimes;
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
        if OCMATAE.freeendtime
            arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
            timehorizon=freepar(OCMATAE.endtimecoordinate);
        else
            arctime=[0 arctime OCMATAE.endtime];
            timehorizon=OCMATAE.endtime;
        end
        if ~OCMATAE.simple
            out.solverinfo.Ycoordinate=OCMATAE.Ycoordinate;
        end
        out.solverinfo.saddlepoint=freepar(OCMATAE.EP.equilibriumcoordinate);
        out.solverinfo.equilibriumcoordinate=OCMATAE.EP.equilibriumcoordinate;
        out.solverinfo.freeparameterindex=OCMATAE.freeparameterindex;
    case 'time'
        arctime=OCMATAE.switchtimes;
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
        %arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        if OCMATAE.freeendtime
            arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
        else
            arctime=[0 arctime OCMATAE.endtime];
        end
        timehorizon=arctime(end);
end
if ~isempty(OCMATAE.freeparameter)
    out.solverinfo.Ycoordinate=OCMATAE.Ycoordinate;
    out.solverinfo.saddlepoint=freepar(OCMATAE.EP.equilibriumcoordinate);
    out.solverinfo.equilibriumcoordinate=OCMATAE.EP.equilibriumcoordinate;
    out.solverinfo.freeparameterindex=OCMATAE.freeparameterindex;
end
out.arcarg=OCMATAE.arcarg;
out.arcinterval=arctime;
out.timehorizon=timehorizon;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.continuationclass='extremal2ep';
out.solverinfo.continuationtype=OCMATCONT.continuationtype;
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.continuationvector=OCMATAE.continuationvector;
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoordinate=OCMATAE.switchtimecoordinate;
out.solverinfo.arcarg=OCMATAE.arcarg;
out.solverinfo.numarc=OCMATAE.numarc;
out.solverinfo.odenum=OCBVP.numode;

out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
if OCMATAE.objectivevaluecalc
    out.solverinfo.objectivevaluecoordinate=OCMATAE.objectivevaluecoordinate;
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

z=[];
modelpar=OCMATAE.parametervalue;
switch OCMATCONT.bvpmethod
    case 'gbvp4c'
        y=zeros(OCBVP.maxnumode,length(tmesh));
        y(OCMATCONT.HE.DDATA.meshvalcoord)=coeff(OCMATCONT.HE.ycoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
end
switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.freeparameterindex)=freemodelparameter;
    otherwise
        if ~isempty(OCMATAE.freeparameter)
            modelpar(OCMATAE.freeparameterindex)=freepar(OCMATAE.continuationcoordinate);
        end

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
global OCMATCONT OCMATAE


if (strcmp(OCMATCONT.continuationtype,'parameter') || ~isempty(OCMATAE.freeparameter)) && ~OCMATAE.simple
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    Y=freepar(OCMATAE.Ycoordinate);
    OCMATCONT.adapted = 1;
    %
    [U,S,V]=svd(OCMATAE.Q0(:,1:OCMATAE.subspacedim)+OCMATAE.Q0(:,OCMATAE.subspacedim+1:end)*Y);
    OCMATAE.Q0=U;
    OCMATAE.Y=zeros(OCMATAE.orthspacedim,OCMATAE.subspacedim);

    freepar(OCMATAE.Ycoordinate)=OCMATAE.Y;

    coeff=[y(:);freepar];
    flag = 1;
elseif (strcmp(OCMATCONT.continuationtype,'parameter') || ~isempty(OCMATAE.freeparameter)) && OCMATAE.simple
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    hatx=freepar(OCMATAE.EP.equilibriumcoordinate);
    Jac=OCMATAE.canonicalsystemjacobian(0,hatx,modelpar,OCMATAE.EP.arcarg);
    if OCMATAE.EP.implicitcontrolnum
        dudx=OCMATAE.dimplicitcontroldx(0,hatx,modelpar,OCMATAE.EP.arcarg);
        Jac=Jac(OCMATAE.statecostatecoordinate,OCMATAE.statecostatecoordinate)+Jac(OCMATAE.statecostatecoordinate,OCMATAE.statecostatecoordinate(end)+(1:OCMATAE.EP.implicitcontrolnum))*dudx;
    end
    OCMATAE.asymptoticmatrix=asymptoticbc(Jac,OCMATAE.pathtype,'c',OCMATCONT.ZeroDeviationTolerance,OCMATCONT.AsymptoticBCMethod);
    flag=1;
else
    flag=0;
end
%-----------------------------------------------------------------
function out=ricatti(Y,J)
global OCMATAE
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=ricatticoefficient(OCMATAE.Q0,J,OCMATAE.subspacedim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);

