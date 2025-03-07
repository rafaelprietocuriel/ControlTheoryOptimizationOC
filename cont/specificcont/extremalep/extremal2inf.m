function out=extremal2inf()

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

switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.continuationparameterindex)=freemodelparameter;
end
transformedtimeshift=0;
if strcmp(OCMATCONT.continuationtype,'time')
    arctime=OCMATAE.switchtimes;
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
    arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
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
    arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
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
if OCMATAE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATAE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(size(J,1),1)];
end
if OCMATAE.exogenousfunction
    J=[[J zeros(size(J,1),OCMATAE.exogenousnumberofstates)];dtds*OCMATAE.exogenousjacobian(t,depvar,modelpar,arcarg)];
end

Jpar=zeros(OCBVP.numode,OCMATCONT.HE.numparameter);
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
        Jpar(1:OCBVP.numode,OCMATAE.switchtimecoordinate(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Jpar(1:OCBVP.numode,OCMATAE.switchtimecoordinate(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        Jpar(1:OCBVP.numode,OCMATAE.switchtimecoordinate(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        if OCMATAE.freeendtime
            Jpar(1:OCBVP.numode,OCMATAE.endtimecoordinate)=dxdt;
            Jpar(1:OCBVP.numode,OCMATAE.switchtimecoordinate(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    end
else
    if OCMATAE.freeendtime
        dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
        Jt=0;
        if OCMATAE.objectivevaluecalc
            dxdt(OCMATAE.objectivevaluecoordinate,:)=OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
            Jt(OCMATAE.objectivevaluecoordinate,:)=OCMATAE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
        end
        if OCMATAE.exogenousfunction
            dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
        end
        Jpar(:,OCMATAE.endtimecoordinate)=dxdt+diffarctime(arc)*(s-arc)*Jt;
    end
end
switch OCMATCONT.continuationtype
    case 'parameter'
        Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
        if isempty(OCMATAE.targetvalue)
            Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATAE.continuationparameterindex);
        else
            Jpar(:,OCMATCONT.HE.numparameter)=OCMATAE.continuationvector.'*Jmodelpar(:,OCMATAE.continuationparameterindex);
        end
    case 'time'
        if ~OCMATAE.autonomous
            Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
        else
            Jt=0;
        end
        if OCMATAE.objectivevaluecalc
            dxdt(OCMATAE.objectivevaluecoordinate,:)=OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
            Jt(OCMATAE.objectivevaluecoordinate,:)=OCMATAE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
        end
        if arc==OCMATAE.continuationtimeindex
            Jpar(1:OCBVP.numode,OCMATAE.continuationcoordinate)=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        elseif any(arc-1==OCMATAE.continuationtimeindex)
            Jpar(1:OCBVP.numode,OCMATAE.continuationcoordinate)=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
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
resconnec=[];
resinit=[];
resend=[];
resricatti=[];
resequilibrium=[];
resasym=[];
userbc=[];

switchtimes=OCMATAE.switchtimes;
switchtimes(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.continuationparameterindex)=freemodelparameter;

        Y=freepar(OCMATAE.Ycoordinate);
        initialstate=OCMATAE.initialpoint(OCMATAE.statecoordinate);
        asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
        resricatti=ricatti(Y,OCMATCONT.monodromymatrix);
        saddlepoint=freepar(OCMATAE.EP.saddlepointcoordinate).';
        resequilibrium=OCMATAE.equilibrium(saddlepoint,modelpar,OCMATAE.EP.arcarg);
        initialstatecoordinate=OCMATAE.statecoordinate;
        endstatecoordinate=[];

    case 'initialstate'
        initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
        asymptoticmatrix=OCMATAE.asymptoticmatrix;
        saddlepoint=OCMATAE.EP.saddlepoint;
        initialstatecoordinate=OCMATAE.depvarcoordinate;
        endvalue=OCMATAE.endpoint(OCMATAE.fixendcoordinate);
        endstatecoordinate=[];
    case 'endstate'
        endstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
        asymptoticmatrix=OCMATAE.asymptoticmatrix;
        saddlepoint=OCMATAE.EP.saddlepoint;
        endstatecoordinate=OCMATAE.depvarcoordinate;
        endvalue=OCMATAE.endpoint(OCMATAE.fixendcoordinate);
        initialstatecoordinate=OCMATAE.fixinitstate;
        initialstate=OCMATAE.initialpoint(OCMATAE.fixinitstate);
    case 'time'
        switchtimes(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        initialstate=OCMATAE.initialpoint(OCMATAE.fixinitstate);
        endvalue=OCMATAE.endpoint(OCMATAE.fixendcoordinate);
        asymptoticmatrix=OCMATAE.asymptoticmatrix;
        saddlepoint=OCMATAE.EP.saddlepoint;
        initialstatecoordinate=OCMATAE.fixinitstate;
        endstatecoordinate=[];
end
if ~isempty(initialstatecoordinate)
    resinit=OCMATAE.bcinitial(depvara(:,OCMATAE.trajectoryindex),initialstatecoordinate,initialstate,modelpar,OCMATAE.EP.arcarg);
end
if ~isempty(endstatecoordinate)
    resend=OCMATAE.bcend(depvarb(:,OCMATAE.trajectoryindex),endstatecoordinate,endstate,modelpar,OCMATAE.EP.arcarg);
end
resinf=OCMATAE.bcinf(depvarb(:,OCMATAE.trajectoryindex(end)),OCMATAE.fixendcoordinate,endvalue,modelpar,OCMATAE.EP.arcarg);
arcarg=OCMATAE.arcarg;

if ~isempty(asymptoticmatrix)
    resasym=OCMATAE.bcasymptotic(depvarb(:,OCMATAE.trajectoryindex(end)),asymptoticmatrix,saddlepoint);
end
if OCMATAE.fixdistance
    resasym=[resasym; ...
        sqrt(sum((saddlepoint(OCMATAE.convergingcoordinate)-depvarb(OCMATAE.convergingcoordinate,OCMATAE.trajectoryindex(end))).^2))-OCMATAE.distance];
end
if OCMATAE.objectivevaluecalc
    resinit=[resinit;depvara(OCMATAE.objectivevaluecoordinate,1)];
end
if OCMATAE.exogenousfunction
    resinit=[resinit; ...
        depvara(OCMATAE.exogenousdynamicscoordinate,1)-OCMATAE.exogenousinitialstates];
end

if OCMATAE.userbc
    userbc=OCMATAE.userfuncbc(depvara,depvarb,modelpar,arcarg,OCMATAE.userbc);
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
end

res=[resinit;resend;resconnec;resasym;resricatti;resequilibrium;resinf;userbc];

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
        arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
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
if 0%~OCMATAE.fixdistance && ~isempty(OCMATAE.convergingcoordinate) &&~strcmp(OCMATCONT.continuationtype,'time')
    switch OCMATCONT.continuationtype
        case 'parameter'
            saddlepoint=freepar(OCMATAE.EP.saddlepointcoordinate).';
        otherwise
            saddlepoint=OCMATAE.EP.saddlepoint;
    end
    yend=y(1:length(saddlepoint),end);

    if isempty(OCMATCONT.OPTIONS.maxdistancecoordinate)
        violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend(OCMATAE.convergingcoordinate)-saddlepoint(OCMATAE.convergingcoordinate))<0;
    else
        violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend(OCMATCONT.OPTIONS.maxdistancecoordinate)-saddlepoint(OCMATCONT.OPTIONS.maxdistancecoordinate))<0;
    end
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

%----------------------------------------------------------------
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
    arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
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
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));

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
                case {'initialstate','parameter','endstate'}
                    out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
                case 'time'
                    if ~isempty(OCMATAE.targetvalue)
                        out=OCMATAE.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
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
[t,y,z,freepar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
switch OCMATCONT.continuationtype
    case {'initialstate','endstate'}
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
        out.solverinfo.Ycoordinate=OCMATAE.Ycoordinate;
    case 'time'
        arctime=OCMATAE.switchtimes;
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.freetimecoordinate);
        if OCMATAE.freeendtime
            arctime=[0 arctime freepar(OCMATAE.endtimecoordinate).'];
        else
            arctime=[0 arctime OCMATAE.endtime];
            arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);

        end
        timehorizon=arctime(end);
end
out.arcarg=OCMATAE.arcarg;
out.arcinterval=arctime;
out.timehorizon=timehorizon;
out.modelparameter=OCMATAE.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremal2inf';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.continuationvector=OCMATAE.continuationvector;
out.solverinfo.pathtype=OCMATAE.pathtype;
%out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoordinate=OCMATAE.switchtimecoordinate;
out.solverinfo.divergingcoordinate=OCMATAE.divergingcoordinate;
out.solverinfo.convergingcoordinate=OCMATAE.convergingcoordinate;
out.solverinfo.fixdistance=OCMATAE.fixdistance;
out.solverinfo.freeendtime=OCMATAE.freeendtime;
out.solverinfo.fixendcoordinate=OCMATAE.fixendcoordinate;

out.solverinfo.arcarg=OCMATAE.arcarg;
out.solverinfo.numarc=OCMATAE.numarc;

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
tmesh=[];
coeff=[];
tangent=[];

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
y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.continuationparameterindex)=freemodelparameter;
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
        save([OCMATAE.basicglobalvarfilename '4extremal2inf'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremal2inf'],'sout','bvpout')
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

%-----------------------------------------------------------------
function out=ricatti(Y,J)
global OCMATAE
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(OCMATAE.Q0,J,OCMATAE.subspacedim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);

