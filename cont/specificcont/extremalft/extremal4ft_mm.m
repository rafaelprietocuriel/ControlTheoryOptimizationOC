function out=extremal4ft_mm()

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
[t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
res=calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)

[t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);

% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,[]},numel(coeff),numJacOpt);
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
arctime=OCMATFTE.initialtimepoints;
arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
end
ct=arctime(OCMATFTE.connectiontimeindex);
part=OCMATFTE.arc2part(arc);
if ~isempty(OCMATFTE.freeparametervector)
    for ii=1:OCMATFTE.numberofstages
        modelpar{ii}(OCMATFTE.freeparameterindex{ii})=freepar(OCMATFTE.freeparametervectorcoordinate{ii});
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofstages
        if ~isempty(OCMATFTE.initialparameter{ii})
            modelpar{ii}(OCMATFTE.continuationindex{ii})=OCMATFTE.initialparameter{ii}+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector{ii};
        end
    end
end

diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));

dtds=diffarctime(arc);
dxdt=dtds*OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
if OCMATFTE.explicitconnectiontime
    dxdt(OCMATFTE.dhamiltoniandctcoord,:)=dtds*OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
end
if OCMATFTE.objectivevaluecalc
    dxdt(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
end
if OCMATFTE.exogenousfunction
    dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=dtds*OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
arctime=OCMATFTE.initialtimepoints;
arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
end
ct=arctime(OCMATFTE.connectiontimeindex);

part=OCMATFTE.arc2part(arc);

diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc); %% test
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc(ones(1,numel(s))));

J=OCMATFTE.canonicalsystemjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
J=dtds*J;
if OCMATFTE.explicitconnectiontime
    J=[J; ...
        dtds*OCMATFTE.dhamiltoniandctjacobian{part}(t,depvar,modelpar{part},arcarg,ct)];
    J=[J zeros(OCMATCONT.DOMAINDDATA{part}(arcindex).numeq+OCMATFTE.connectiontimenumber-1,OCMATFTE.connectiontimenumber)];
end
if OCMATFTE.objectivevaluecalc
    if OCMATFTE.explicitconnectiontime
        J=[J; ...
            [dtds*OCMATFTE.objectivefunctionjacobian{part}(t,depvar,modelpar{part},arcarg,ct) zeros(1,OCMATFTE.connectiontimenumber)]];
        J=[J zeros(OCMATCONT.DOMAINDDATA{part}(arcindex).numeq+OCMATFTE.connectiontimenumber,1)];
    else
        J=[J; ...
            dtds*OCMATFTE.objectivefunctionjacobian{part}(t,depvar,modelpar{part},arcarg,ct)];
        J=[J zeros(OCMATCONT.DOMAINDDATA{part}(arcindex).numeq,1)];
    end
end
if OCMATFTE.exogenousfunction
    J=[J; ...
        [dtds*OCMATFTE.exogenousjacobian{part}(t,depvar,modelpar{part},arcarg,ct) zeros(OCMATFTE.exogenousnumberofstates,OCMATFTE.connectiontimenumber+1)]];
    J=[J zeros(OCMATCONT.DOMAINDDATA{part}(arcindex).numeq+OCMATFTE.connectiontimenumber+OCMATFTE.exogenousnumberofstates,OCMATFTE.exogenousnumberofstates)];
end

if OCMATFTE.explicitconnectiontime
    Jpar=zeros(OCMATCONT.DOMAINDDATA{part}(arcindex).numeq+OCMATFTE.connectiontimenumber,OCMATCONT.HE.numparameter);
else
    Jpar=zeros(OCMATCONT.DOMAINDDATA{part}(arcindex).numeq,OCMATCONT.HE.numparameter);
end
if OCMATFTE.exogenousfunction
    Jpar=[Jpar;zeros(OCMATFTE.exogenousnumberofstates,OCMATCONT.HE.numparameter)];
end
if OCMATFTE.timepointscharacterization(arc+1)==1
    dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.explicitconnectiontime
        dxdt(OCMATFTE.dhamiltoniandctcoord,:)=OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.exogenousfunction
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    Jpar(:,OCMATFTE.timepointsfree4parcoordinate(arc+1==OCMATFTE.timepointsfreeindex))=dxdt;
end
if OCMATFTE.timepointscharacterization(arc)==1
    dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.explicitconnectiontime
        dxdt(OCMATFTE.dhamiltoniandctcoord,:)=OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.exogenousfunction
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    Jpar(:,OCMATFTE.timepointsfree4parcoordinate(arc==OCMATFTE.timepointsfreeindex))=-dxdt;
end
if OCMATFTE.timepointscharacterization(arc+1)==-1
    dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
        Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.explicitconnectiontime
        dxdt(OCMATFTE.dhamiltoniandctcoord,:)=OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.exogenousfunction
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    Jpar(:,OCMATFTE.timepointsoptimize4parcoordinate(arc+1==OCMATFTE.timepointsoptimizedindex))=dxdt;%+Jt;
end
if OCMATFTE.timepointscharacterization(arc)==-1
    dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
        Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.explicitconnectiontime
        dxdt(OCMATFTE.dhamiltoniandctcoord,:)=OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.exogenousfunction
        dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    Jpar(:,OCMATFTE.timepointsoptimize4parcoordinate(arc==OCMATFTE.timepointsoptimizedindex))=-dxdt;%+Jt;
end
if OCMATFTE.continuationtype==1
    if  OCMATFTE.timepointscharacterization(arc)==2
        dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
        if OCMATFTE.objectivevaluecalc
            dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
            Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{part}(t,depvar,modelpar{part},arcarg,ct);
        end
        if OCMATFTE.explicitconnectiontime
            dxdt(OCMATFTE.dhamiltoniandctcoord,:)=OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
        end
        if OCMATFTE.exogenousfunction
            dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
        end
        idx=find(OCMATFTE.continuationindex==arc);
        Jpar(:,OCMATFTE.continuationcoordinate)=Jpar(:,OCMATFTE.continuationcoordinate)-OCMATFTE.continuationvector(idx)*dxdt;
    end
    if  OCMATFTE.timepointscharacterization(arc+1)==2
        dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
        if OCMATFTE.objectivevaluecalc
            dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
            Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{part}(t,depvar,modelpar{part},arcarg,ct);
        end
        if OCMATFTE.explicitconnectiontime
            dxdt(OCMATFTE.dhamiltoniandctcoord,:)=OCMATFTE.dhamiltoniandct{part}(t,depvar,modelpar{part},arcarg,ct);
        end
        if OCMATFTE.exogenousfunction
            dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics{part}(t,depvar,modelpar{part},arcarg,ct);
        end
        idx=find(OCMATFTE.continuationindex==arc+1);
        Jpar(:,OCMATFTE.continuationcoordinate)=Jpar(:,OCMATFTE.continuationcoordinate)+OCMATFTE.continuationvector(idx)*dxdt;
    end
    %     if arc==OCMATCONT.HE.numarc  && any(OCMATFTE.continuationindex==OCMATCONT.HE.numarc+1)
    %         dxdt=OCMATFTE.canonicalsystem{part}(t,depvar,modelpar{part},arcarg,ct);
    %         if OCMATFTE.objectivevaluecalc
    %             dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{part}(t,depvar,modelpar{part},arcarg,ct);
    %             Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{part}(t,depvar,modelpar{part},arcarg,ct);
    %         end
    %         idx=find(OCMATFTE.continuationindex==OCMATCONT.HE.numarc+1);
    %         Jpar(:,OCMATFTE.continuationcoordinate)=Jpar(:,OCMATFTE.continuationcoordinate)-OCMATFTE.continuationvector(idx)*dxdt;
    %     end
end
if OCMATFTE.freeparameter
    Jmodelpar=dtds*OCMATFTE.canonicalsystemparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
    Jpar(OCMATFTE.statecostatecoord{part},OCMATFTE.parametercoord)=Jmodelpar(:,OCMATFTE.parameterindex{part});
    if OCMATFTE.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
        Jpar(OCMATFTE.objectivevaluecoord,OCMATFTE.parametercoord)=Jobjmodelpar(:,OCMATFTE.parameterindex{part});
    end
    if OCMATFTE.explicitconnectiontime
        Jdhdctmodelpar=dtds*OCMATFTE.dhamiltoniandctparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
        Jpar(OCMATFTE.dhamiltoniandctcoord,OCMATFTE.parametercoord)=Jdhdctmodelpar(:,OCMATFTE.parameterindex{part});
    end
    if OCMATFTE.exogenousfunction
        Jexfmodelpar=dtds*OCMATFTE.exogenousparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
        Jpar(OCMATFTE.exogenousdynamicscoordinate,OCMATFTE.parametercoord)=Jexfmodelpar(:,OCMATFTE.parameterindex{part});
    end
end
if OCMATFTE.continuationtype==2
    Jmodelpar=dtds*OCMATFTE.canonicalsystemparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
    if OCMATFTE.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.explicitconnectiontime
        Jdhdctmodelpar=dtds*OCMATFTE.dhamiltoniandctparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.exogenousfunction
        Jexfmodelpar=dtds*OCMATFTE.exogenousparameterjacobian{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if length(OCMATFTE.continuationvector{part})==1
        Jpar(OCMATFTE.statecostatecoord{part},end)=OCMATFTE.continuationvector{part}*Jmodelpar(:,OCMATFTE.continuationindex{part});
        if OCMATFTE.objectivevaluecalc
            Jpar(OCMATFTE.objectivevaluecoord,end)=OCMATFTE.continuationvector{part}*Jobjmodelpar(:,OCMATFTE.continuationindex{part});
        end
        if OCMATFTE.explicitconnectiontime
            Jpar(OCMATFTE.dhamiltoniandctcoord,end)=OCMATFTE.continuationvector{part}*Jdhdctmodelpar(:,OCMATFTE.continuationindex{part});
        end
        if OCMATFTE.exogenousfunction
            Jpar(OCMATFTE.exogenousdynamicscoordinate,end)=OCMATFTE.continuationvector{part}*Jexfmodelpar(:,OCMATFTE.continuationindex{part});
        end
    else
        for ii=1:length(OCMATFTE.continuationvector{part})
            Jpar(OCMATFTE.statecostatecoord{part},end)=Jpar(OCMATFTE.statecostatecoord{part},end)+OCMATFTE.continuationvector{part}(ii)*Jmodelpar(:,OCMATFTE.continuationindex{part}(ii));
            if OCMATFTE.objectivevaluecalc
                Jpar(OCMATFTE.objectivevaluecoord,end)=Jpar(OCMATFTE.objectivevaluecoord,end)+OCMATFTE.continuationvector{part}(ii)*Jobjmodelpar(:,OCMATFTE.continuationindex{part}(ii));
            end
            if OCMATFTE.explicitconnectiontime
                Jpar(OCMATFTE.dhamiltoniandctcoord,end)=Jpar(OCMATFTE.dhamiltoniandctcoord,end)+OCMATFTE.continuationvector{part}(ii)*Jdhdctmodelpar(:,OCMATFTE.continuationindex{part}(ii));
            end
            if OCMATFTE.exogenousfunction
                Jpar(OCMATFTE.exogenousdynamicscoordinate,end)=Jpar(OCMATFTE.exogenousdynamicscoordinate,end)+OCMATFTE.continuationvector{part}(ii)*Jexfmodelpar(:,OCMATFTE.continuationindex{part}(ii));
            end
        end
    end
end
if OCMATFTE.explicitconnectiontime
    d2xdtdc=dtds*OCMATFTE.derivativeconnectiontime{part}(t,depvar,modelpar{part},arcarg,ct);
    if OCMATFTE.objectivevaluecalc
        d2xdtdc(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunctionderivativeconnectiontime{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    if OCMATFTE.exogenousfunction
        d2xdtdc(OCMATFTE.exogenousdynamicscoordinate,:)=dtds*OCMATFTE.exogenousdct{part}(t,depvar,modelpar{part},arcarg,ct);
    end
    d2xdtdc(OCMATFTE.dhamiltoniandctcoord,:)=dtds*OCMATFTE.d2hamiltoniandct2{part}(t,depvar,modelpar{part},arcarg,ct);
    idx=find(OCMATFTE.optimalswitchingtime(1:end-1));
    Jpar(:,OCMATFTE.timepointsoptimize4parcoordinate)=Jpar(:,OCMATFTE.timepointsoptimize4parcoordinate)+d2xdtdc(:,idx);
end
%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP

resconnec=[];
resinit=[];
restrans=[];
restarget=[];
respartconnect=[];

arctime=OCMATFTE.initialtimepoints;
arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
end
ct=arctime(OCMATFTE.connectiontimeindex);

if ~isempty(OCMATFTE.freeparametervector)
    for ii=1:OCMATFTE.numberofstages
        modelpar{ii}(OCMATFTE.freeparameterindex{ii})=freepar(OCMATFTE.freeparametervectorcoordinate{ii});
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofstages
        if ~isempty(OCMATFTE.initialparameter{ii})
            modelpar{ii}(OCMATFTE.continuationindex{ii})=OCMATFTE.initialparameter{ii}+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector{ii};
        end
    end
end


if OCMATFTE.objectivevaluecalc
    OVal=0;
end
if OCMATFTE.explicitconnectiontime
    dHdctval=zeros(OCMATFTE.connectiontimenumber,1);
end
for arc=1:OCMATCONT.HE.numarc
    part=OCMATFTE.arc2part(arc);
    arcarg=OCMATCONT.HE.arcarg;
    if arc==1
        initialstate=OCMATFTE.initialstate;
        if OCMATFTE.continuationtype==0
            initialstate=initialstate+freepar(end)*OCMATFTE.continuationvector;
        end
        restarget=depvara(OCMATFTE.statecoordinate{part},1)-initialstate;
        if OCMATFTE.userbc
            restarget=[restarget; ...
                OCMATFTE.userfunctionbc{arc}(depvara,depvarb,modelpar,ct,OCMATCONT.HE.arcarg,arc,part)];
        end

    end
    if OCMATFTE.objectivevaluecalc && arc==OCMATCONT.HE.numarc
        OVal=OVal+OCMATFTE.salvagevalue{part}(arctime(arc+1),depvarb(OCMATFTE.statecostatecoord{part},arc),modelpar{part},OCMATCONT.HE.arcarg(arc),ct);
    end
    if OCMATFTE.explicitconnectiontime && any(arc==OCMATFTE.partstructure)
        dHdctval=dHdctval+OCMATFTE.dsalvagedct{part}(arctime(arc+1),depvarb(OCMATFTE.statecostatecoord{part},arc),modelpar{part},OCMATCONT.HE.arcarg(arc),ct);
    end
    if OCMATFTE.timepointscharacterization(arc+1)==1
        resconnec=[resconnec; ...
            OCMATFTE.reset{part}(depvara(OCMATFTE.statecostatecoord{part},:),depvarb(OCMATFTE.statecostatecoord{part},:),modelpar{part},arctime,OCMATCONT.HE.arcarg,OCMATFTE.edge,arc,ct); ...
            OCMATFTE.guard{part}(depvara(OCMATFTE.statecostatecoord{part},:),depvarb(OCMATFTE.statecostatecoord{part},:),modelpar{part},arctime,OCMATCONT.HE.arcarg,OCMATFTE.edge,arc,ct)];
        if OCMATFTE.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.objectivevaluecoord,arc)-depvara(OCMATFTE.objectivevaluecoord,arc+1)];
        end
        if OCMATFTE.exogenousfunction
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.exogenousdynamicscoordinate,arc)-depvara(OCMATFTE.exogenousdynamicscoordinate,arc+1)];
        end
    end
    if part<OCMATFTE.numberofstages && any(arc==OCMATFTE.partstructure) % arc at the end of a part
        respartconnect=[respartconnect; ...
            OCMATFTE.bcconnectingparts{part}(depvara(OCMATFTE.statecostatecoord{part},:),depvarb(OCMATFTE.statecostatecoord{part},:),modelpar{part},ct,arcarg,arc)];
        if  OCMATFTE.timepointscharacterization(arc+1)==-1
            respartconnect=[respartconnect; ...
                OCMATFTE.bcoptimalconnectingparts{part}(depvara,depvarb,modelpar,ct,OCMATCONT.HE.arcarg,arc,part)];
        end
        if OCMATFTE.objectivevaluecalc
            respartconnect=[respartconnect; ...
                depvara(OCMATFTE.objectivevaluecoord,arc+1)-depvarb(OCMATFTE.objectivevaluecoord,arc)];
        end
        if OCMATFTE.explicitconnectiontime
            respartconnect=[respartconnect; ...
                depvara(OCMATFTE.dhamiltoniandctcoord,arc+1)-depvarb(OCMATFTE.dhamiltoniandctcoord,arc)];
        end
        if OCMATFTE.exogenousfunction
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.exogenousdynamicscoordinate,arc)-depvara(OCMATFTE.exogenousdynamicscoordinate,arc+1)];
        end
    end
    if arc==OCMATCONT.HE.numarc
        if OCMATFTE.timepointscharacterization(arc+1)==-1
            % end time optimally chosen
        end
        restrans=[restrans; ...
            OCMATFTE.bctransversality{part}(arctime(end),depvarb(OCMATFTE.statecostatecoord{part},arc),modelpar{part},arcarg,ct)];
        if OCMATFTE.objectivevaluecalc
            resinit=[resinit; ...
                depvara(OCMATFTE.objectivevaluecoord,1)-OVal];
        end
        if OCMATFTE.explicitconnectiontime
            resinit=[resinit; ...
                depvara(OCMATFTE.dhamiltoniandctcoord,1)-dHdctval];
        end
        if OCMATFTE.exogenousfunction
            resinit=[resinit; ...
                depvara(OCMATFTE.exogenousdynamicscoordinate,1)-OCMATFTE.exogenousinitialstates];
        end
    elseif arc==1 && OCMATFTE.timepointscharacterization(arc)==-1
        % initial time optimally chosen
    end

end

res=[resinit;restarget;restrans;resconnec;respartconnect];

%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];
%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATFTE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    part=OCMATFTE.arc2part(arc);
    arctime=OCMATFTE.initialtimepoints;
    arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
    arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
    if OCMATFTE.continuationtype==1
        arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
    end
    if ~isempty(OCMATFTE.freeparametervector)
        for ii=1:OCMATFTE.numberofstages
            modelpar{ii}(OCMATFTE.freeparameterindex{ii})=freepar(OCMATFTE.freeparametervectorcoordinate{ii});
        end
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATFTE.testadmissibility{part}(t,sol.y(idx,:),modelpar{part},arcarg,ct);
    else
        eqcoord=domainddata{part}(arcindex).eqcoord;
        [constr labelS]=OCMATFTE.testadmissibility{part}(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar{part},arcarg,ct);
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
        %infoS(counter).constraintvalue=diffarctime(arc);
        infoS(counter).minval=min(diffarctime(:));
        b=min([b infoS(counter).minval]);
    end
end


%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATFTE OCMATCONT
[s,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
arctime=OCMATFTE.initialtimepoints;
arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
end
diffarctime=diff(arctime);
for arc=1:numel(OCMATCONT.HE.arcarg)
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(arc)-diffarctime(arc)*(arc-1));
end
h=OCMATFTE.plotcontinuation{1}(t,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE
idx=[];
if isempty(coeff)
    return
end
if OCMATFTE.findoptimalswitchingtime
    [t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
    arctime=OCMATFTE.initialtimepoints;
    arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
    arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
    if OCMATFTE.continuationtype==1
        arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
    end
    arc=OCMATFTE.findoptimalswitchingtime-1;
    part=OCMATFTE.arc2part(arc);
    depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
    depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
    arcarg=OCMATCONT.HE.arcarg;
    out=OCMATFTE.bcoptimalconnectingparts{part}(depvara,depvarb,modelpar,arctime,arcarg,arc,part);
    fprintf(1,' Difference Hamiltonian: %g\n',out);
else
    out=coeff(OCMATCONT.HE.contparametercoord);
    fprintf(1,' Continuation parameter: %g\n',out);
end
if OCMATFTE.findoptimalparameter
    tpar=tangent(end);
    tangent=tangent/tpar;
    tangent=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
    fprintf(1,' Derivative value      : %g\n',tangent(OCMATFTE.objectivevaluecoord,end));
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

            if OCMATFTE.findoptimalswitchingtime
                [t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
                arctime=OCMATFTE.initialtimepoints;
                arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
                arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
                if OCMATFTE.continuationtype==1
                    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
                end
                ct=arctime(OCMATFTE.connectiontimeindex);
                arc=OCMATFTE.findoptimalswitchingtime-1;
                part=OCMATFTE.arc2part(arc);
                depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                arcarg=OCMATCONT.HE.arcarg;
                out=OCMATFTE.bcoptimalconnectingparts{part}(depvara,depvarb,modelpar,ct,arcarg,arc,part);
            elseif OCMATFTE.findoptimalparameter
                tpar=tangent(end);
                tangent=tangent/tpar;
                tangent=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
                out=tangent(OCMATFTE.objectivevaluecoord,end);
                OCMATFTE.derivative=out;
            elseif OCMATFTE.hitvalue
                [t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);
                arctime=OCMATFTE.initialtimepoints;
                arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
                arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
                if OCMATFTE.continuationtype==1
                    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
                end
                ct=arctime(OCMATFTE.connectiontimeindex);
                %arc=OCMATFTE.findoptimalswitchingtime-1;
                %part=OCMATFTE.arc2part(arc);
                depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                arcarg=OCMATCONT.HE.arcarg;
                out=OCMATFTE.hitvaluefunc(depvara,depvarb,modelpar,arcarg,ct);
            else
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
            end
        otherwise
            out=[];
            failed=1;

    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE OCBVP
[t,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
if ~isempty(OCMATFTE.freeparametervector)
    for ii=1:OCMATFTE.numberofstages
        modelpar{ii}(OCMATFTE.freeparameterindex{ii})=freepar(OCMATFTE.freeparametervectorcoordinate{ii});
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofstages
        if ~isempty(OCMATFTE.initialparameter{ii})
            modelpar{ii}(OCMATFTE.continuationindex{ii})=OCMATFTE.initialparameter{ii}+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector{ii};
        end
    end
end

arctime=OCMATFTE.initialtimepoints;
arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
end
out.arcinterval=arctime;
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=arctime(1);
out.timehorizon=arctime(end);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremal4ft_mm';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.continuationindex=OCMATFTE.continuationindex;
out.solverinfo.partstructure=OCMATFTE.partstructure;
out.solverinfo.partposition=OCMATFTE.partposition;
out.solverinfo.optimalswitchingtime=OCMATFTE.optimalswitchingtime;
out.solverinfo.connectingtimepointsindex=OCMATFTE.connectiontimeindex;
out.solverinfo.continuationtype=OCMATFTE.continuationtype;
out.solverinfo.numberofstages=OCMATFTE.numberofstages;
out.solverinfo.dhamiltoniandctcoord=OCMATFTE.dhamiltoniandctcoord;
out.solverinfo.objectivevaluecoord=OCMATFTE.objectivevaluecoord;
out.solverinfo.exogenousdynamicscoordinate=OCMATFTE.exogenousdynamicscoordinate;
out.solverinfo.userbc=OCMATFTE.userbc;

if OCMATFTE.continuationtype==0
elseif OCMATFTE.continuationtype==1
    out.solverinfo.continuationvector=OCMATFTE.continuationvector;
    out.solverinfo.initialtime4continuation=OCMATFTE.initialtime4continuation;
elseif OCMATFTE.continuationtype==2
    out.solverinfo.continuationindex=OCMATFTE.continuationindex;
    out.solverinfo.initialparameter=OCMATFTE.initialparameter;
    out.solverinfo.continuationvector=OCMATFTE.continuationvector;
end

switch OCMATCONT.bvpmethod
    case 'bvp4c'
        out.solverinfo.yp=out.yp;
        out=rmfield(out,{'yp'});
end

% out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
% out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
% if ~isempty(OCMATFTE.fixendstatecoord)
%     out.solverinfo.fixendstatecoord=OCMATFTE.fixendstatecoord;
% end

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
function [tmesh,y,z,freepar,modelpar,ct]=drearr(tmesh,coeff)
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
if ~isempty(OCMATFTE.freeparametervector)
    for ii=1:OCMATFTE.numberofstages
        modelpar{ii}(OCMATFTE.freeparameterindex{ii})=freepar(OCMATFTE.freeparametervectorcoordinate{ii});
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofstages
        if ~isempty(OCMATFTE.initialparameter{ii})
            modelpar{ii}(OCMATFTE.continuationindex{ii})=OCMATFTE.initialparameter{ii}+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.continuationvector{ii};
        end
    end
end
arctime=OCMATFTE.initialtimepoints;
arctime(OCMATFTE.timepointsoptimizedindex)=freepar(OCMATFTE.timepointsoptimize4parcoordinate);
arctime(OCMATFTE.timepointsfreeindex)=freepar(OCMATFTE.timepointsfree4parcoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationindex)=OCMATFTE.initialtime4continuation+freepar(end)*OCMATFTE.continuationvector;
end
ct=arctime(OCMATFTE.connectiontimeindex);


%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATFTE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATFTE=OCMATFTE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATFTE.basicglobalvarfilename '4extremal4ft_mm'],'MODELINFO')
    end
    save([OCMATFTE.basicresultfilename '4extremal4ft_mm'],'sout','bvpout')
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


