function out=extremalp4ft()

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

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
modelpar(OCMATFTE.continuationindex)=OCMATFTE.initialparametervalue+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector;
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
end
if ~OCMATFTE.optimalhorizon
    arctime=OCMATFTE.initarcinterval;
    arctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.freeswitchingtimecoordinate);
    %arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATFTE.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.vfreetimecoord);
    vdiffarctime=diff(varctime);
    vdts=vdiffarctime(arc);
    dxdt(OCMATFTE.variationaldynamicscoordinate,:)=dtds*OCMATFTE.variationaldynamics(t,depvar,modelpar,arcarg)+vdts*OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
end
if OCMATFTE.objectivevaluecalc
    dxdt(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
    if OCMATFTE.includevariationalobjectivevalue
        dxdt(OCMATFTE.variationalobjectivevaluecoord,:)=dtds*OCMATFTE.variationalobjectivefunction(t,depvar,modelpar,arcarg)+vdts*OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
    end
end
if OCMATFTE.exogenousfunction
    dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=dtds*OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
    if OCMATFTE.variationalcalculation
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=dxdt(OCMATFTE.exogenousdynamicscoordinate,:)+vdts*OCMATFTE.exogenousvariationaldynamics(t,depvar,modelpar,arcarg);
    end
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
modelpar(OCMATFTE.continuationindex)=OCMATFTE.initialparametervalue+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector;
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
end
if ~OCMATFTE.optimalhorizon
    arctime=OCMATFTE.initarcinterval;
    arctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.freeswitchingtimecoordinate);
    %arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
diffarctime=diff(arctime);
%arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATFTE.JX;
if OCMATFTE.variationalcalculation
    J0=OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg); % J0 is reused
    J(OCMATFTE.dFDXcoord1,OCMATFTE.dFDXcoord2)=dtds*J0;
else
    J(OCMATFTE.dFDXcoord1,OCMATFTE.dFDXcoord2)=dtds*OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
end
if OCMATFTE.objectivevaluecalc
    J(OCMATFTE.dFODXcoord1(1),OCMATFTE.dFODXcoord2)=dtds*OCMATFTE.objectivefunctionjacobian(t,depvar,modelpar,arcarg);
    if OCMATFTE.includevariationalobjectivevalue
        J(OCMATFTE.dFODXcoord1(2),OCMATFTE.dFODXcoord2)=dtds*OCMATFTE.variationalobjectivefunction(t,depvar,modelpar,arcarg);
    end
end
if OCMATFTE.exogenousfunction
    J(OCMATFTE.dFEDXcoord1,OCMATFTE.dFEDXcoord2)=dtds*OCMATFTE.exogenousjacobian(t,depvar,modelpar,arcarg);
end
if OCMATFTE.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.vfreetimecoord);
    vdiffarctime=diff(varctime);
    vdts=vdiffarctime(arc);
    %J(OCMATFTE.dFVDXcoord1,OCMATFTE.ODEcoord)=dtds*[OCMATFTE.variationaljacobian(t,depvar,modelpar,arcarg) OCMATFTE.dFDE]+vdts*[J0 OCMATFTE.dFVDX OCMATFTE.dFDO OCMATFTE.dFDE];
    J(OCMATFTE.dFVDXcoord1,OCMATFTE.ODEcoord)=dtds*OCMATFTE.variationaljacobian(t,depvar,modelpar,arcarg)+vdts*[J0 OCMATFTE.dFVDX OCMATFTE.dFDO OCMATFTE.dFDE];
    if OCMATFTE.exogenousfunction
        J(OCMATFTE.dFEDXcoord1,OCMATFTE.ODEcoord)=J(OCMATFTE.dFEDXcoord1,OCMATFTE.ODEcoord)+[dtds*OCMATFTE.exogenousjacobian4variationalargument(t,depvar,modelpar,arcarg) OCMATFTE.dFEDO OCMATFTE.dFEDE];%+[vdts*OCMATFTE.exogenousvariationaldynamicsjacobian(t,depvar,modelpar,arcarg) OCMATFTE.dFEDO OCMATFTE.dFEDE];
    end
end
Jpar=OCMATFTE.Jpar;
Jmodelpar=OCMATFTE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATFTE.autonomous
        Jt=OCMATFTE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATFTE.variationalcalculation
        vdxdt(OCMATFTE.variationaldynamicscoordinate,:)=dxdt;
        dxdt(OCMATFTE.variationaldynamicscoordinate,:)=OCMATFTE.variationaldynamics(t,depvar,modelpar,arcarg);
        Jt(OCMATFTE.variationaldynamicscoordinate,:)=0;%OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
        if OCMATFTE.variationalcalculation
            if OCMATFTE.includevariationalobjectivevalue
                dxdt(OCMATFTE.variationalobjectivevaluecoord,:)=OCMATFTE.variationalobjectivefunction(t,depvar,modelpar,arcarg);
                vdxdt(OCMATFTE.objectivevaluecoord,:)=0;
                Jt(OCMATFTE.variationalobjectivevaluecoord,:)=OCMATFTE.variationalobjectivefunctionderivativetime(t,depvar,modelpar,arcarg);
                vdxdt(OCMATFTE.variationalobjectivevaluecoord,:)=0;
            else
                vdxdt(OCMATFTE.objectivevaluecoord,:)=0;
            end
        end
    end
    if OCMATFTE.exogenousfunction
        if OCMATFTE.variationalcalculation
            vdxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousvariationaldynamics(t,depvar,modelpar,arcarg);
        end
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
        Jt(OCMATFTE.exogenousdynamicscoordinate,:)=0;
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATFTE.ODEcoord,OCMATFTE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if OCMATFTE.variationalcalculation
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.vfreetimecoord(arc))=vdxdt;
        end
        if arc>1
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
            if OCMATFTE.variationalcalculation
                Jpar(OCMATFTE.ODEcoord,OCMATFTE.vfreetimecoord(arc-1))=-vdxdt;
            end
        end
    else
        Jpar(OCMATFTE.ODEcoord,OCMATFTE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        if OCMATFTE.optimalhorizon
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        end
        if OCMATFTE.variationalcalculation
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.vfreetimecoord(arc-1))=-vdxdt;
        end
    end
else
    if ~OCMATFTE.autonomous
        Jt=OCMATFTE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATFTE.optimalhorizon
        dxdt=OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(OCMATFTE.ODEcoord,OCMATFTE.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc)*Jt;
    end
end

if length(OCMATFTE.continuationindex)==1
    Jpar(OCMATFTE.statecostatecoord,end)=OCMATFTE.continuationvector*dtds*Jmodelpar(:,OCMATFTE.continuationindex);
else
    Jpar(OCMATFTE.statecostatecoord,end)=sum(repmat(OCMATFTE.continuationvector,size(Jmodelpar(:,OCMATFTE.continuationindex),1),1).*dtds.*Jmodelpar(:,OCMATFTE.continuationindex),2);
end
if OCMATFTE.variationalcalculation
    Jvarmodelpar=dtds*OCMATFTE.variationalparameterjacobian(t,depvar,modelpar,arcarg);
    if length(OCMATFTE.continuationindex)==1
        Jpar(OCMATFTE.variationaldynamicscoordinate,end)=OCMATFTE.continuationvector*(Jvarmodelpar(:,OCMATFTE.continuationindex)+vdts*Jmodelpar(:,OCMATFTE.continuationindex));
    else
        %Jpar(OCMATFTE.variationaldynamicscoordinate,end)=sum(repmat(OCMATFTE.continuationvector,size(Jvarmodelpar(:,OCMATFTE.continuationindex),1),1).*(Jvarmodelpar(:,OCMATFTE.continuationindex)+vdts*Jmodelpar(:,OCMATFTE.continuationindex)),2);
    end
end
if OCMATFTE.objectivevaluecalc
    Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian(t,depvar,modelpar,arcarg);
    if OCMATFTE.includevariationalobjectivevalue
        Jvarobjmodelpar=dtds*OCMATFTE.variationalobjectivefunctionparameterjacobian(t,depvar,modelpar,arcarg);
    end
    if length(OCMATFTE.continuationindex)==1
        Jpar(OCMATFTE.objectivevaluecoord,end)=OCMATFTE.continuationvector*Jobjmodelpar(:,OCMATFTE.continuationindex);
        if OCMATFTE.includevariationalobjectivevalue
            Jpar(OCMATFTE.variationalobjectivevaluecoord,end)=OCMATFTE.continuationvector*Jvarobjmodelpar(:,OCMATFTE.continuationindex);
        end
    else
        Jpar(OCMATFTE.objectivevaluecoord,end)=sum(repmat(OCMATFTE.continuationvector,size(Jobjmodelpar(:,OCMATFTE.continuationindex),1),1).*Jobjmodelpar(:,OCMATFTE.continuationindex),2);
    end
end
if OCMATFTE.exogenousfunction
    Jexfmodelpar=dtds*OCMATFTE.exogenousparameterjacobian(t,depvar,modelpar,arcarg);
    if length(OCMATFTE.continuationindex)==1
        Jpar(OCMATFTE.exogenousdynamicscoordinate,end)=OCMATFTE.continuationvector*Jexfmodelpar(:,OCMATFTE.continuationindex);
    else
        Jpar(OCMATFTE.exogenousdynamicscoordinate,end)=sum(repmat(OCMATFTE.continuationvector,size(Jexfmodelpar(:,OCMATFTE.continuationindex),1),1).*Jexfmodelpar(:,OCMATFTE.continuationindex),2);
    end
end
if ~isempty(OCMATFTE.freeparameter)
    Jpar(OCMATFTE.statecostatecoord,OCMATFTE.freeparametercoordinate)=dtds*Jmodelpar(:,OCMATFTE.freeparameter);
    if OCMATFTE.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian(t,depvar,modelpar,arcarg);
        Jpar(OCMATFTE.objectivevaluecoord,OCMATFTE.freeparametercoordinate)=Jobjmodelpar(:,OCMATFTE.freeparameter);
    end
    if OCMATFTE.exogenousfunction
        Jexfmodelpar=dtds*OCMATFTE.exogenousparameterjacobian(t,depvar,modelpar,arcarg);
        Jpar(OCMATFTE.exogenousdynamicscoordinate,OCMATFTE.freeparametercoordinate)=Jexfmodelpar(:,OCMATFTE.freeparameter);
    end
    if OCMATFTE.variationalcalculation
        Jpar(OCMATFTE.variationaldynamicscoordinate,OCMATFTE.freeparametercoordinate)=Jvarmodelpar(:,OCMATFTE.freeparameter)+vdts*Jmodelpar(:,OCMATFTE.freeparameter);%OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
end


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
modelpar(OCMATFTE.continuationindex)=OCMATFTE.initialparametervalue+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector;
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
end
arctime=OCMATFTE.initarcinterval;
arctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.freeswitchingtimecoordinate);
switchtimes=arctime(2:end-1);
if ~OCMATFTE.optimalhorizon
    timehorizon=OCMATFTE.truncationtime;
else
    timehorizon=freepar(OCMATFTE.optimalhorizoncoord);
    if OCMATFTE.variationalcalculation
        vtimehorizon=freepar(OCMATFTE.voptimalhorizoncoord);
    end
end
resconnec=[];
resinit=[];
resobj=[];
userbc=[];
if ~isempty(OCMATFTE.fixinitstatecoord)
    resinit=OCMATFTE.bcinitial(depvara,OCMATFTE.fixinitstatecoord,OCMATFTE.initstate,modelpar,OCMATCONT.HE.arcarg(1));
    if ~isempty(OCMATFTE.freeinitstatecoord)
        initstate=OCMATFTE.initstate;
        initstate(OCMATFTE.freeinitstatecoord)=freepar(OCMATFTE.freeinitstateidx);
        resinit=[resinit; ...
            OCMATFTE.bcinitial(depvara,OCMATFTE.freeinitstatecoord,initstate,modelpar,OCMATCONT.HE.arcarg(1))];
    end
    if OCMATFTE.variationalcalculation
        resinit=[resinit; ...
            OCMATFTE.variationalbcinitial(depvara,[],[],modelpar,OCMATCONT.HE.arcarg(1));];
    end
end
if OCMATFTE.constraintidx
    val=OCMATFTE.testadmissibility(OCMATFTE.initialtime,depvara(:,1),modelpar,OCMATCONT.HE.arcarg(1));
    resinit=[resinit; ...
        val(OCMATFTE.constraintidx)];
end
if OCMATFTE.userbc
    userbc=OCMATFTE.userbcfunc([OCMATFTE.initialtime switchtimes(:).' timehorizon],depvara,depvarb,modelpar,OCMATCONT.HE.arcarg);
end
if OCMATFTE.stateconstraint
    jumparg=freepar(OCMATFTE.entrytimecoordinate);
end
if OCMATFTE.variationalcalculation
    vswitchtimes=freepar(OCMATFTE.vfreetimecoord);
    if OCMATFTE.stateconstraint
        vjumparg=freepar(OCMATFTE.ventrytimecoordinate);
    end
end

if OCMATFTE.transversalityconditioncs
    restrans=OCMATFTE.bctransversalitysc(timehorizon,depvarb,jumparg(end),modelpar,OCMATFTE.jumpid(end));
    if OCMATFTE.variationalcalculation
        restrans=[restrans; ...
            OCMATFTE.variationalbctransversalitysc(timehorizon,depvarb,[jumparg(end) vjumparg(end)],modelpar,OCMATFTE.jumpid(end));];
    end
else
    restrans=OCMATFTE.bctransversality(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end));
    if OCMATFTE.variationalcalculation
        restrans=[restrans; ...
            OCMATFTE.variationalbctransversality(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end));];
    end
end
if ~isempty(OCMATFTE.fixendstatecoord)
    restrans(OCMATFTE.fixendstatecoord,1)=depvarb(OCMATFTE.fixendstatecoord,end)-OCMATFTE.endstate;
end
if OCMATFTE.optimalhorizon
    restrans=[restrans; ...
        OCMATFTE.bcoptimalhorizon(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end))];
    if OCMATFTE.variationalcalculation
    end
end
if ~isempty(OCMATFTE.followobjectivevalue)
    resobj=OCMATFTE.followobjectivevalue-depvarb(OCMATFTE.objectivevaluecoord,end);
end

if OCMATFTE.objectivevaluecalc
    OVal=OCMATFTE.salvagevalue(timehorizon,depvarb(:,end),modelpar,OCMATCONT.HE.arcarg(end));
    resinit=[resinit;depvara(OCMATFTE.objectivevaluecoord,1)-OVal];
    if OCMATFTE.includevariationalobjectivevalue
        varOVal=OCMATFTE.variationalsalvagevalue(timehorizon,depvarb(:,end),modelpar,OCMATCONT.HE.arcarg(end));
        resinit=[resinit;depvara(OCMATFTE.variationalobjectivevaluecoord,1)-varOVal];
    end
end
if OCMATFTE.exogenousfunction
    resinit=[resinit; ...
        depvara(OCMATFTE.exogenousdynamicscoordinate,1)-OCMATFTE.exogenousinitialstates(timehorizon,depvara,depvarb,modelpar,OCMATCONT.HE.arcarg)];
%        depvara(OCMATFTE.exogenousdynamicscoordinate,1)-OCMATFTE.exogenousinitialstates];
end

for ii=1:numel(OCMATCONT.HE.arcarg)-1
    if OCMATFTE.stateconstraint && OCMATFTE.entryindex(ii)
        if OCMATFTE.switchingtimetype(ii+1)
            resconnec=[resconnec; ...
                OCMATFTE.bcstateconstraint(depvara,depvarb,modelpar,jumparg(OCMATFTE.entryindex(ii)),switchtimes,OCMATFTE.jumpid(ii),OCMATCONT.HE.edge,ii); ...
                OCMATFTE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
        else
            resconnec=[resconnec;
                OCMATFTE.bcstateconstraint(depvara,depvarb,modelpar,jumparg(OCMATFTE.entryindex(ii)),switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
        end
    else
        if OCMATFTE.switchingtimetype(ii+1)
            resconnec=[resconnec;
                OCMATFTE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
                OCMATFTE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
        else
            resconnec=[resconnec;
                OCMATFTE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
        end
    end
    if OCMATFTE.variationalcalculation
        if OCMATFTE.stateconstraint && OCMATFTE.entryindex(ii)
            if OCMATFTE.switchingtimetype(ii+1)
                resconnec=[resconnec; ...
                    OCMATFTE.variationalbcstateconstraint(depvara,depvarb,modelpar,[jumparg(OCMATFTE.entryindex(ii)) vjumparg(OCMATFTE.entryindex(ii))],switchtimes,vswitchtimes,OCMATFTE.jumpid(ii),OCMATCONT.HE.edge,ii); ...
                    OCMATFTE.variationalguard(depvara,depvarb,modelpar,switchtimes,vswitchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
            else
                resconnec=[resconnec;
                    OCMATFTE.variationalbcstateconstraint(depvara,depvarb,modelpar,[jumparg(OCMATFTE.entryindex(ii)) vjumparg(OCMATFTE.entryindex(ii))],switchtimes,vswitchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
            end
        else
            if OCMATFTE.switchingtimetype(ii+1)
                resconnec=[resconnec;
                    OCMATFTE.variationalreset(depvara,depvarb,modelpar,switchtimes,vswitchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
                    OCMATFTE.variationalguard(depvara,depvarb,modelpar,switchtimes,vswitchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
            else
                resconnec=[resconnec;
                    OCMATFTE.variationalreset(depvara,depvarb,modelpar,switchtimes,vswitchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
            end
        end
    end
    if OCMATFTE.objectivevaluecalc
        resconnec=[resconnec; ...
            depvarb(OCMATFTE.objectivevaluecoord,ii)-depvara(OCMATFTE.objectivevaluecoord,ii+1)];
        if OCMATFTE.includevariationalobjectivevalue
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.variationalobjectivevaluecoord,ii)-depvara(OCMATFTE.variationalobjectivevaluecoord,ii+1)];
        end
    end
    if OCMATFTE.exogenousfunction
        resconnec=[resconnec; ...
            depvarb(OCMATFTE.exogenousdynamicscoordinate,ii)-depvara(OCMATFTE.exogenousdynamicscoordinate,ii+1)];
    end
end

res=[resinit;resconnec;restrans;resobj;userbc];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATFTE OCBVP

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if ~OCMATFTE.optimalhorizon
        arctime=OCMATFTE.initarcinterval;     
        arctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.freeswitchingtimecoordinate);     
        %arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime]; 
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
        [constr,labelS]=OCMATFTE.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        %eqcoord=domainddata(arcindex).eqcoord;
        [constr,labelS]=OCMATFTE.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    if OCMATFTE.stateconstraint
        jumparg=freepar(OCMATFTE.entrytimecoordinate);
        violationmat=jumparg<-OCMATCONT.OPTIONS.admissibletol;
        if any(violationmat)
            counter=counter+1;
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='stateconstraint';
            infoS(counter).cols=1;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=jumparg;
            infoS(counter).minval=min(jumparg);
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
    arctime=OCMATFTE.initarcinterval;     
    arctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.freeswitchingtimecoordinate);     
    %arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime]; 
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
diffarctime=diff(arctime);
for arc=1:numel(OCMATCONT.HE.arcarg)
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(arc)-diffarctime(arc)*(arc-1));
end
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATFTE.plotcontinuation(t,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
%drawnow
%figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
if isempty(OCMATFTE.hitvalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
else
    fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
    fprintf(1,' Difference to hitvalue: %g\n',OCMATFTE.hitvaluefunc(t,y,modelpar,OCMATCONT.HE.arcarg,OCMATFTE.hitvalue));
end
if OCMATFTE.findoptimalparameter
    %tpar=tangent(end);
    %tangent=tangent/tpar;
    %tangent=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
    fprintf(1,' Derivative value      : %g\n',y(OCMATFTE.variationalobjectivevaluecoord,end));
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
[t,y,z,freepar]=drearr(tmesh,coeff);

failed=[];
for ii=id
    switch ii
        case 1
            if OCMATFTE.findoptimalparameter
                out=y(OCMATFTE.variationalobjectivevaluecoord,end);
                OCMATFTE.derivative=out;
            elseif isempty(OCMATFTE.hitvalue)
                out=1-coeff(end);
            else
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                out=OCMATFTE.hitvaluefunc(t,y,modelpar,OCMATCONT.HE.arcarg,OCMATFTE.hitvalue);
            end

        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
if ~OCMATFTE.optimalhorizon
    arctime=OCMATFTE.initarcinterval;     
    arctime(OCMATFTE.freeswitchingtimeindex)=freepar(OCMATFTE.freeswitchingtimecoordinate);     
    %arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime]; 
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
out.arcinterval=arctime;
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
out.solverinfo.conttype='extremalp4ft';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATFTE.switchtimecoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.solverinfo.continuationindex=OCMATFTE.continuationindex;
out.solverinfo.continuationcoordinate=OCMATFTE.continuationcoordinate;
out.solverinfo.continuationvector=OCMATFTE.continuationvector;
out.solverinfo.initialparametervaluevalue=OCMATFTE.initialparametervalue;
out.solverinfo.fixendstatecoord=OCMATFTE.fixendstatecoord;
out.solverinfo.fixinitstatecoord=OCMATFTE.fixinitstatecoord;
out.solverinfo.freeparametercoordinate=OCMATFTE.freeparametercoordinate;
out.solverinfo.freeswitchingtimecoordinate=OCMATFTE.freeswitchingtimecoordinate;
out.solverinfo.constraintidx=OCMATFTE.constraintidx;
out.solverinfo.objectivevaluecalc=OCMATFTE.objectivevaluecalc;
out.solverinfo.objectivevaluecoord=OCMATFTE.objectivevaluecoord;
out.solverinfo.optimalhorizon=OCMATFTE.optimalhorizon;
out.solverinfo.followobjectivevalue=OCMATFTE.followobjectivevalue;
out.solverinfo.stateconstraint=OCMATFTE.stateconstraint;
if OCMATFTE.variationalcalculation
    out.solverinfo.vfreetimecoord=OCMATFTE.vfreetimecoord;
end
if ~isempty(OCMATFTE.hitvalue)
    out.solverinfo.hitvalue=OCMATFTE.hitvaluefunc(t,y,modelpar,OCMATCONT.HE.arcarg,OCMATFTE.hitvalue);
end
if OCMATFTE.stateconstraint
    out.solverinfo.entrytimecoordinate=OCMATFTE.entrytimecoordinate;
    out.solverinfo.entryindex=OCMATFTE.entryindex;
    out.solverinfo.jumpid=OCMATFTE.jumpid;
    if OCMATFTE.variationalcalculation
        out.solverinfo.ventrytimecoordinate=OCMATFTE.ventrytimecoordinate;
    end
end
if OCMATFTE.exogenousfunction
    out.solverinfo.exogenousdynamicscoordinate=OCMATFTE.exogenousdynamicscoordinate;
end

if OCMATFTE.optimalhorizon
    out.solverinfo.optimalhorizoncoord=OCMATFTE.optimalhorizoncoord;
end

if OCMATFTE.findoptimalparameter
    out.solverinfo.dodpar=y(OCMATFTE.variationalobjectivevaluecoord,end);
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
global OCMATCONT OCMATFTE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATFTE.PD_phi=p'/norm(p);
        OCMATFTE.PD_psi=Q(:,end);
        s.data.phi=OCMATFTE.PD_phi(:);
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
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATFTE OCBVP

modelpar=OCMATFTE.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        % handle case of pure state constraints
    otherwise
end
modelpar(OCMATFTE.continuationindex)=OCMATFTE.initialparametervalue+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector;
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
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
        save([OCMATFTE.basicglobalvarfilename '4extremalp4ft'],'MODELINFO')
    end
    save([OCMATFTE.basicresultfilename '4extremalp4ft'],'sout','bvpout')
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