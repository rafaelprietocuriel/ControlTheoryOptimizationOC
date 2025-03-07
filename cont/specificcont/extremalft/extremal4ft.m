function out=extremal4ft()

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

function [F,J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
if OCMATFTE.freeparameter
    modelpar(OCMATFTE.parameterindex)=freepar(OCMATFTE.parametercoord);
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATFTE.variationalcalculation
    dxdt(OCMATFTE.variationaldynamicscoordinate,:)=dtds*OCMATFTE.variationaldynamics(t,depvar,modelpar,arcarg);
end
if OCMATFTE.objectivevaluecalc
    dxdt(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
end
if OCMATFTE.exogenousfunction
    dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=dtds*OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
end
if OCMATFTE.exogenousfunction
    dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=dtds*OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
% if OCMATFTE.freeparameter
%     modelpar(OCMATFTE.parameterindex{ii})=freepar(OCMATFTE.parametercoord);
% end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATFTE.objectivevaluecalc
    J=[J;dtds*OCMATFTE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
end
if OCMATFTE.exogenousfunction
    J=[J;dtds*OCMATFTE.exogenousjacobian(t,depvar,modelpar,arcarg)];
end
J=[J,OCMATFTE.Jext];
Jpar=OCMATFTE.Jpar;
if OCMATCONT.HE.numarc>1
    dxdt=OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATFTE.autonomous
        Jt=OCMATFTE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATFTE.exogenousfunction
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
        Jt(OCMATFTE.exogenousdynamicscoordinate,:)=0;
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATFTE.ODEcoord,OCMATFTE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        if arc>1
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
    else
        Jpar(OCMATFTE.ODEcoord,OCMATFTE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        if OCMATFTE.optimalhorizon
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc)*Jt;
        end
    end
else
    if OCMATFTE.optimalhorizon
        dxdt=OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
        if OCMATFTE.objectivevaluecalc
            dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
            if ~OCMATFTE.autonomous
                Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
            else
                Jt=0;
            end
        end
        if OCMATFTE.exogenousfunction
            dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
            Jt(OCMATFTE.exogenousdynamicscoordinate,:)=0;
        end
        Jpar(OCMATFTE.ODEcoord,OCMATFTE.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc)*Jt;
    end
end
if OCMATFTE.freeparameter
    Jmodelpar=dtds*OCMATFTE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
    Jpar(OCMATFTE.statecostatecoord,OCMATFTE.parametercoord)=Jmodelpar(:,OCMATFTE.parameterindex);
    if OCMATFTE.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian(t,depvar,modelpar,arcarg);
        Jpar(OCMATFTE.objectivevaluecoord,OCMATFTE.parametercoord)=Jobjmodelpar(:,OCMATFTE.parameterindex);
    end
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP

switchtimes=freepar(OCMATFTE.switchtimecoord);
if ~OCMATFTE.optimalhorizon
    timehorizon=OCMATFTE.truncationtime;
else
    timehorizon=freepar(OCMATFTE.optimalhorizoncoord);
end
if OCMATFTE.freeparameter
    modelpar(OCMATFTE.parameterindex)=freepar(OCMATFTE.parametercoord);
end
resconnec=[];
resend=[];
resuser=[];
if OCMATFTE.continitstate
    initialstate=OCMATFTE.startvalue+freepar(end)*OCMATFTE.continuationvector;
    resinit=OCMATFTE.bcinitial(depvara,OCMATFTE.targetcoordinate,initialstate,modelpar,OCMATCONT.HE.arcarg(1));
end
if ~isempty(OCMATFTE.freeinitstatecoord)
    resinit=[resinit; ...
        OCMATFTE.bcfreeinitial(depvara,OCMATFTE.freeinitstatecoord,modelpar,OCMATCONT.HE.arcarg(1),freepar(OCMATFTE.initstateinequalitycoord))];
end
if OCMATFTE.stateconstraint
    jumparg=freepar(OCMATFTE.entrytimecoordinate);
end

if OCMATFTE.transversalityconditioncs
    restrans=OCMATFTE.bctransversalitysc(timehorizon,depvarb,jumparg(end),modelpar,OCMATFTE.jumpid(end));
else
    restrans=OCMATFTE.bctransversality(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end));
end
if ~isempty(OCMATFTE.fixendstatecoord)
    restrans(OCMATFTE.fixendstatecoord,1)=depvarb(OCMATFTE.fixendstatecoord,end)-OCMATFTE.endstate;
end
if ~OCMATFTE.continitstate
    restrans(OCMATFTE.targetcoordinate,1)=depvarb(OCMATFTE.targetcoordinate,end)-endstate;
end
if OCMATFTE.variationalcalculation
    resinit=[resinit; ...
        depvara(OCMATFTE.variationaldynamicscoordinate,1)-OCMATFTE.variationalinitialstates];
end

if OCMATFTE.objectivevaluecalc
    OVal=OCMATFTE.salvagevalue(timehorizon,depvarb(:,end),modelpar,OCMATCONT.HE.arcarg(end));
    resinit=[resinit;depvara(OCMATFTE.objectivevaluecoord,1)-OVal];
end
if OCMATFTE.exogenousfunction
    resinit=[resinit; ...
        depvara(OCMATFTE.exogenousdynamicscoordinate,1)-OCMATFTE.exogenousinitialstates];
end
if OCMATFTE.optimalhorizon
    restrans=[restrans; ...
        OCMATFTE.bcoptimalhorizon(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end))];
end
if ~isempty(OCMATFTE.userbc)
    resuser=[resuser OCMATFTE.userbc(depvara,depvarb,modelpar,OCMATCONT.HE.arcarg,OCMATFTE.userbcvalue);];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    if OCMATFTE.stateconstraint && OCMATFTE.entryindex(ii)
        resconnec=[resconnec; ...
            OCMATFTE.bcstateconstraint(depvara,depvarb,modelpar,jumparg(OCMATFTE.entryindex(ii)),switchtimes,OCMATFTE.jumpid(ii),OCMATCONT.HE.edge,ii); ...
            OCMATFTE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    else
        resconnec=[resconnec;
            OCMATFTE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATFTE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
    if OCMATFTE.variationalcalculation
        resconnec=[resconnec; ...
            depvarb(OCMATFTE.variationaldynamicscoordinate,ii)-depvara(OCMATFTE.variationaldynamicscoordinate,ii+1)];
    end
    if OCMATFTE.objectivevaluecalc
        resconnec=[resconnec; ...
            depvarb(OCMATFTE.objectivevaluecoord,ii)-depvara(OCMATFTE.objectivevaluecoord,ii+1)];
    end
    if OCMATFTE.exogenousfunction
        resconnec=[resconnec; ...
            depvarb(OCMATFTE.exogenousdynamicscoordinate,ii)-depvara(OCMATFTE.exogenousdynamicscoordinate,ii+1)];
    end
end
res=[resinit;resend;resconnec;restrans;resuser];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
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
switchtimes=freepar(OCMATFTE.switchtimecoord);
if OCMATFTE.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATFTE.jumpcostateindex)=freepar(OCMATFTE.jumpcostatecoord);
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
[Japart Jbpart Jpartmp]=OCMATFTE.bcjacobianinitial(depvara,freepar,OCMATCONT.HE.arcarg(1),OCMATFTE.targetcoordinate,OCMATFTE.continuationvector);
if OCMATFTE.objectivevaluecalc
    Japart(OCMATFTE.objectivevaluecoord,OCMATFTE.objectivevaluecoord)=1;
end

Jpar(1:OCMATCONT.HE.numinitialcondition,:)=Jpartmp;
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numinitialcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;
for arc=1:OCMATCONT.HE.numarc-1
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    Jaresetpart=OCMATFTE.jacobianreset(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,-1);
    Jbresetpart=OCMATFTE.jacobianreset(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,0);
    Jaguardpart=OCMATFTE.jacobianguard(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,-1);
    Jbguardpart=OCMATFTE.jacobianguard(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,0);

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

[Japart Jbpart  Jpartmp]=OCMATFTE.bcjacobianasymptotic(depvarb,OCMATCONT.HE.arcarg(OCMATCONT.HE.numarc),OCMATFTE.asymptoticmatrix',OCMATFTE.saddlepoint);
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numendcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;


%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATFTE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if ~OCMATFTE.optimalhorizon
        arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
    else
        arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATFTE.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        %eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATFTE.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    if isfinite(OCMATFTE.maxhorizon)
        violationmat=OCMATFTE.maxhorizon-arctime(end)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
        if violationmat
            counter=counter+1;
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='maxhorizon';
            infoS(counter).cols=1;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=OCMATFTE.maxhorizon-arctime(end);
            infoS(counter).minval=OCMATFTE.maxhorizon-arctime(end);
            b=min([b infoS(counter).minval]);
        end
    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATFTE OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
diffarctime=diff(arctime);
for arc=1:numel(OCMATCONT.HE.arcarg)
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(arc)-diffarctime(arc)*(arc-1));
end
h=OCMATFTE.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
% drawnow
% figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE
idx=[];
if isempty(coeff)
    return
end
if isempty(OCMATFTE.hitvalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
else
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    if ~OCMATFTE.optimalhorizon
        arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
    else
        arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
    end
    fprintf(1,' distance hit function: %g\n',OCMATFTE.hitvalue-OCMATFTE.hitvaluefunc(arctime,y,modelpar,OCMATCONT.HE.arcarg));
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATFTE
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
global OCMATCONT OCMATFTE
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
global OCMATCONT OCMATFTE

failed=[];
for ii=id
    switch ii
        case 1
            if isempty(OCMATFTE.hitvalue)
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
            else
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                if ~OCMATFTE.optimalhorizon
                    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
                else
                    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
                end

                out=OCMATFTE.hitvalue-OCMATFTE.hitvaluefunc(arctime,y,modelpar,OCMATCONT.HE.arcarg);
            end

        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
if OCMATFTE.freeparameter
    modelpar(OCMATFTE.parameterindex)=freepar(OCMATFTE.parametercoord);
end
out.arcinterval=arctime;
%out.arcinterval=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATFTE.initialtime;
out.timehorizon=arctime(end);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremal4ft';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATFTE.switchtimecoord;
out.solverinfo.optimalhorizoncoord=OCMATFTE.optimalhorizoncoord;
out.solverinfo.stateconstraint=OCMATFTE.stateconstraint;
if OCMATFTE.stateconstraint
    out.solverinfo.entrytimecoordinate=OCMATFTE.entrytimecoordinate;    
    out.solverinfo.entryindex=OCMATFTE.entryindex;    
    out.solverinfo.jumpid=OCMATFTE.jumpid;    
end

out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
if ~isempty(OCMATFTE.fixendstatecoord)
    out.solverinfo.fixendstatecoord=OCMATFTE.fixendstatecoord;
end
if OCMATFTE.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATFTE.objectivevaluecoord;
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
if OCMATFTE.stateconstraint
%    out.solverinfo.jumpcostatecoord=OCMATFTE.jumpcostatecoord;
end

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATFTE
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
        s.data.laecoefficient=calchnf_LP(tmesh,coeff,tangent,@operatoreq,q,p,J);
        s.data.sol=formatsolution(tmesh,coeff,tangent);
        s.msg =sprintf('Limit extremal');
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
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATFTE OCBVP

modelpar=OCMATFTE.parametervalue;

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
if OCMATFTE.freeparameter
    modelpar(OCMATFTE.parameterindex)=freepar(OCMATFTE.parametercoord);
end


%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATFTE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATFTE=OCMATFTE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATFTE.basicglobalvarfilename '4extremal4ft'],'MODELINFO')
    end
    save([OCMATFTE.basicresultfilename '4extremal4ft'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATFTE

pathname=OCMATFTE.datapath();

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATFTE

discretizationdata=OCMATFTE.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATFTE
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
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
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end