function out=indifferencesolution4ft()

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
global OCMATINDIF OCMATCONT OCBVP
solutionindex=OCMATINDIF.solutionindex(arc);
modelpar=modelpar{solutionindex};
if OCMATINDIF.continuationtype==2
    modelpar(OCMATINDIF.continuationindex)=OCMATINDIF.initialparameter+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
end
if OCMATINDIF.freeparametervector
    if isempty(OCMATINDIF.freeinitialparametervalue)
        modelpar(OCMATINDIF.freeparameterindex)=freepar(OCMATINDIF.freeparametervectorcoordinate);
    else
        modelpar(OCMATINDIF.freeparameterindex)=OCMATINDIF.freeinitialparametervalue+freepar(OCMATINDIF.freeparametervectorcoordinate)*OCMATINDIF.freeparameterdirection;
    end
end
if OCMATINDIF.static
    modelpar(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{solutionindex});
end
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if OCMATINDIF.continuationtype==1
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.initialendtime+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
end
if ~isempty(OCMATINDIF.freetimevector)
    arctime(end)=freepar(OCMATINDIF.freetimevectorcoordinate);
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATINDIF.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(2:end-1)=freepar(OCMATINDIF.vfreetimecoord{solutionindex});
    vdiffarctime=diff(varctime);
    vdts=vdiffarctime(relarc);
    dxdt(OCMATINDIF.variationaldynamicscoordinate,:)=dtds*OCMATINDIF.variationaldynamics(t,depvar,modelpar,arcarg)+vdts*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
end
if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
    if OCMATINDIF.variationalcalculation
        dxdt(OCMATINDIF.variationalobjectivevaluecoord,:)=dtds*OCMATINDIF.variationalobjectivefunction(t,depvar,modelpar,arcarg);
    end
end
if OCMATINDIF.exogenousfunction
    dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=dtds*OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
    if OCMATINDIF.variationalcalculation
        dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)+vdts*OCMATINDIF.exogenousvariationaldynamics(t,depvar,modelpar,arcarg);
    end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
solutionindex=OCMATINDIF.solutionindex(arc);
modelpar=modelpar{solutionindex};
if OCMATINDIF.static
    modelpar(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{solutionindex});
end
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if OCMATINDIF.continuationtype==1
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.initialendtime+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
end
if ~isempty(OCMATINDIF.freetimevector)
    arctime(end)=freepar(OCMATINDIF.freetimevectorcoordinate);
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
%arcindex=OCMATCONT.HE.arcindex(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATINDIF.JX;
if OCMATINDIF.variationalcalculation
    J0=OCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg); % J0 is reused
    J(OCMATINDIF.dFDXcoord1,OCMATINDIF.dFDXcoord2)=dtds*J0;
else
    J(OCMATINDIF.dFDXcoord1,OCMATINDIF.dFDXcoord2)=dtds*OCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
end
if OCMATINDIF.objectivevaluecalc
    J(OCMATINDIF.dFODXcoord1,OCMATINDIF.dFODXcoord2)=dtds*OCMATINDIF.objectivefunctionjacobian(t,depvar,modelpar,arcarg);
end
if OCMATINDIF.exogenousfunction
    J(OCMATINDIF.dFEDXcoord1,OCMATINDIF.dFEDXcoord2)=dtds*OCMATINDIF.exogenousjacobian(t,depvar,modelpar,arcarg);
end
if OCMATINDIF.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(2:end-1)=freepar(OCMATINDIF.vfreetimecoord{solutionindex});
    vdiffarctime=diff(varctime);
    vdts=vdiffarctime(relarc);
    J(OCMATINDIF.dFVDXcoord1,OCMATINDIF.ODEcoord)=dtds*OCMATINDIF.variationaljacobian(t,depvar,modelpar,arcarg)+vdts*[J0 OCMATINDIF.dFVDX OCMATINDIF.dFDO OCMATINDIF.dFDE];
    if OCMATINDIF.exogenousfunction
        J(OCMATINDIF.dFEDXcoord1,OCMATINDIF.ODEcoord)=J(OCMATINDIF.dFEDXcoord1,OCMATINDIF.ODEcoord)+[OCMATINDIF.dFEDX dtds*OCMATINDIF.exogenousjacobian4variationalargument(t,depvar,modelpar,arcarg) OCMATINDIF.dFEDO OCMATINDIF.dFEDE]+[vdts*OCMATINDIF.exogenousvariationaldynamicsjacobian(t,depvar,modelpar,arcarg) OCMATINDIF.dFEDO OCMATINDIF.dFEDE];
    end
end
Jpar=OCMATINDIF.Jpar;
if OCMATINDIF.numarc(solutionindex)>1
    dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATINDIF.autonomous
        Jt=OCMATINDIF.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATINDIF.variationalcalculation
        vdxdt(OCMATINDIF.variationaldynamicscoordinate,:)=dxdt;
        dxdt(OCMATINDIF.variationaldynamicscoordinate,:)=OCMATINDIF.variationaldynamics(t,depvar,modelpar,arcarg);
        Jt(OCMATINDIF.variationaldynamicscoordinate,:)=0;%OCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATINDIF.exogenousfunction
        if OCMATINDIF.variationalcalculation
            vdxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousvariationaldynamics(t,depvar,modelpar,arcarg);
        end
        dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
        Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
    end
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-relarc)*Jt;
        if OCMATINDIF.variationalcalculation
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.vfreetimecoord{solutionindex}(relarc))=vdxdt;
        end
        if relarc>1
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc-1)*Jt);
            if OCMATINDIF.variationalcalculation
                Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.vfreetimecoord{solutionindex}(relarc-1))=-vdxdt;
            end
        end
    else
        Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc-1)*Jt);
        if OCMATINDIF.variationalcalculation
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.vfreetimecoord{solutionindex}(relarc-1))=-vdxdt;
        end
        if OCMATINDIF.continuationtype==1
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            if OCMATINDIF.objectivevaluecalc
                dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
                %Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
            end
            if OCMATINDIF.exogenousfunction
                dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
                %Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
            end
            Jpar(:,OCMATINDIF.continuationcoordinate)=OCMATINDIF.continuationvector*dxdt;
            if OCMATINDIF.variationalcalculation
                Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.continuationcoordinate)=OCMATINDIF.continuationvector*vdxdt;
            end
        end
    end
else
    if OCMATINDIF.continuationtype==1
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        if OCMATINDIF.objectivevaluecalc
            dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
            %Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
        end
        if OCMATINDIF.exogenousfunction
            dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
            Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
        end
        Jpar(:,OCMATINDIF.continuationcoordinate)=OCMATINDIF.continuationvector*dxdt;
    end
end

if OCMATINDIF.continuationtype==2
    Jmodelpar=OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
    if length(OCMATINDIF.continuationindex)==1
        Jpar(OCMATINDIF.statecostatecoord,end)=OCMATINDIF.continuationvector*dtds*Jmodelpar(:,OCMATINDIF.continuationindex);
    else
        Jpar(OCMATINDIF.statecostatecoord,end)=sum(repmat(OCMATINDIF.continuationvector,size(Jmodelpar(:,OCMATINDIF.continuationindex),1),1).*dtds.*Jmodelpar(:,OCMATINDIF.continuationindex),2);
    end
    if OCMATINDIF.variationalcalculation
        Jvarmodelpar=dtds*OCMATINDIF.variationalparameterjacobian(t,depvar,modelpar,arcarg);
        if length(OCMATINDIF.continuationindex)==1
            Jpar(OCMATINDIF.variationaldynamicscoordinate,end)=OCMATINDIF.continuationvector*(Jvarmodelpar(:,OCMATINDIF.continuationindex)+vdts*Jmodelpar(:,OCMATINDIF.continuationindex));
        else
            %Jpar(OCMATINDIF.variationaldynamicscoordinate,end)=sum(repmat(OCMATINDIF.continuationvector,size(Jvarmodelpar(:,OCMATINDIF.continuationindex),1),1).*(Jvarmodelpar(:,OCMATINDIF.continuationindex)+vdts*Jmodelpar(:,OCMATINDIF.continuationindex)),2);
        end
    end
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian(t,depvar,modelpar,arcarg);
        if length(OCMATINDIF.continuationindex)==1
            Jpar(OCMATINDIF.objectivevaluecoord,end)=OCMATINDIF.continuationvector*Jobjmodelpar(:,OCMATINDIF.continuationindex);
        else
            Jpar(OCMATINDIF.objectivevaluecoord,end)=sum(repmat(OCMATINDIF.continuationvector,size(Jobjmodelpar(:,OCMATINDIF.continuationindex),1),1).*Jobjmodelpar(:,OCMATINDIF.continuationindex),2);
        end
    end
    if OCMATINDIF.exogenousfunction
        Jexfmodelpar=dtds*OCMATINDIF.exogenousparameterjacobian(t,depvar,modelpar,arcarg);
        if length(OCMATINDIF.continuationindex)==1
            Jpar(OCMATINDIF.exogenousdynamicscoordinate,end)=OCMATINDIF.continuationvector*Jexfmodelpar(:,OCMATINDIF.continuationindex);
        else
            Jpar(OCMATINDIF.exogenousdynamicscoordinate,end)=sum(repmat(OCMATINDIF.continuationvector,size(Jexfmodelpar(:,OCMATINDIF.continuationindex),1),1).*Jobjmodelpar(:,OCMATINDIF.continuationindex),2);
        end
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    Jmodelpar=OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
    if isempty(OCMATINDIF.freeinitialparametervalue)
        freeJmodelpar=dtds*Jmodelpar(:,OCMATINDIF.freeparameterindex);
    else
        freeJmodelpar=sum(repmat(OCMATINDIF.freeparameterdirection,OCMATINDIF.statecostatecoord(end),1).*dtds*Jmodelpar(:,OCMATINDIF.freeparameterindex),2);
    end
    Jpar(OCMATINDIF.statecostatecoord,OCMATINDIF.freeparametervectorcoordinate)=freeJmodelpar;
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian(t,depvar,modelpar,arcarg);
        if isempty(OCMATINDIF.freeinitialparametervalue)
            freeJobjmodelpar=Jobjmodelpar(:,OCMATINDIF.freeparameterindex);
        else
            freeJobjmodelpar=sum(OCMATINDIF.freeparameterdirection.*Jobjmodelpar(:,OCMATINDIF.freeparameterindex),2);
        end
        Jpar(OCMATINDIF.objectivevaluecoord,OCMATINDIF.freeparametervectorcoordinate)=freeJobjmodelpar;
    end
    if OCMATINDIF.exogenousfunction
        Jexfmodelpar=dtds*OCMATINDIF.exogenousparameterjacobian(t,depvar,modelpar,arcarg);
        if isempty(OCMATINDIF.freeinitialparametervalue)
            freeJexjmodelpar=Jexfmodelpar(:,OCMATINDIF.freeparameterindex);
        else
            freeJexjmodelpar=sum(OCMATINDIF.freeparameterdirection.*Jexfmodelpar(:,OCMATINDIF.freeparameterindex),2);
        end
        Jpar(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.freeparametervectorcoordinate)=freeJexjmodelpar;
    end
    if OCMATINDIF.variationalcalculation
        if isempty(OCMATINDIF.freeinitialparametervalue)
            freeJvajmodelpar=Jvarmodelpar(:,OCMATINDIF.freeparameterindex)+vdts*Jmodelpar(:,OCMATINDIF.freeparameterindex);
        else
            freeJvajmodelpar=sum(OCMATINDIF.freeparameterdirection.*Jvarmodelpar(:,OCMATINDIF.freeparameterindex)+vdts*Jmodelpar(:,OCMATINDIF.freeparameterindex),2);
        end
        Jpar(OCMATINDIF.variationaldynamicscoordinate,OCMATINDIF.freeparametervectorcoordinate)=freeJvajmodelpar;
    end
end
if ~isempty(OCMATINDIF.variationalcalculation) && OCMATINDIF.static
    Jpar(OCMATINDIF.statecostatecoord,OCMATINDIF.staticparametercoordinate{solutionindex})=dtds*Jmodelpar(:,OCMATINDIF.staticparameterindex);
    Jpar(OCMATINDIF.variationaldynamicscoordinate,OCMATINDIF.staticparametercoordinate{solutionindex})=Jvarmodelpar(:,OCMATINDIF.staticparameterindex)+vdts*Jmodelpar(:,OCMATINDIF.staticparameterindex);
end


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT
resconnec=[];
residpt=[];
resinit=[];
restrans=[];
restarget=[];
userbc=[];
modelpar0=modelpar;
% if OCMATINDIF.continuationtype==2
%     modelpar(OCMATINDIF.continuationindex)=OCMATINDIF.initialparameter+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
% end
% if OCMATINDIF.freeparametervector
%     if isempty(OCMATINDIF.freeinitialparametervalue)
%         modelpar(OCMATINDIF.freeparameterindex)=freepar(OCMATINDIF.freeparametervectorcoordinate);
%     else
%         modelpar(OCMATINDIF.freeparameterindex)=OCMATINDIF.freeinitialparametervalue+freepar(OCMATINDIF.freeparametervectorcoordinate)*OCMATINDIF.freeparameterdirection;
%     end
% end

O=zeros(1,OCMATINDIF.indifferenceorder);
if OCMATINDIF.continuationtype==1
    endtime=OCMATINDIF.initialendtime+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
else
    endtime=OCMATINDIF.endtime;
end
if ~isempty(OCMATINDIF.freetimevector)
    endtime=freepar(OCMATINDIF.freetimevectorcoordinate);
end

for ii=1:OCMATINDIF.indifferenceorder
    modelpar=modelpar0{ii};
    if OCMATINDIF.continuationtype==2
        modelpar(OCMATINDIF.continuationindex)=OCMATINDIF.initialparameter+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
    end
    if OCMATINDIF.freeparametervector
        if isempty(OCMATINDIF.freeinitialparametervalue)
            modelpar(OCMATINDIF.freeparameterindex)=freepar(OCMATINDIF.freeparametervectorcoordinate);
        else
            modelpar(OCMATINDIF.freeparameterindex)=OCMATINDIF.freeinitialparametervalue+freepar(OCMATINDIF.freeparametervectorcoordinate)*OCMATINDIF.freeparameterdirection;
        end
    end
    if OCMATINDIF.static
        modelpar(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{ii});
    end
    switchtimes=freepar(OCMATINDIF.switchtimecoord{ii});
    arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii});
    if ii==1
        initialstate=OCMATINDIF.startvalue;
        if ~isempty(OCMATINDIF.freestatevector)
            initialstate=initialstate+sum(repmat(freepar(OCMATINDIF.freestatevectorcoordinate),1,OCMATINDIF.statecostatecoordinate(end)).*OCMATINDIF.freestatevector,2);
        end
        if OCMATINDIF.continuationtype==0
            initialstate=initialstate+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
        end
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate,1)-initialstate];
    end
    if OCMATINDIF.variationalcalculation
        resinit=[resinit; ...
            OCMATINDIF.variationalbcinitial(depvara(:,OCMATINDIF.arccoord{ii}(1)),[],[],modelpar,arcarg(1))];
    end
    if ii<OCMATINDIF.indifferenceorder
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(ii+1))-depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(ii))];
    end
    if OCMATINDIF.objectivevaluecalc
        OVal=OCMATINDIF.salvagevalue(endtime,depvarb(OCMATINDIF.statecostatecoord,OCMATINDIF.arccoord{ii}(end)),modelpar,arcarg(end));
        resinit=[resinit; ...
            depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(1))-OVal];
        if OCMATINDIF.variationalcalculation
            varOVal=OCMATINDIF.variationalsalvagevalue(endtime,depvarb(:,OCMATINDIF.arccoord{ii}(end)),modelpar,OCMATCONT.HE.arcarg(end));
            resinit=[resinit;depvara(OCMATINDIF.variationalobjectivevaluecoord,OCMATINDIF.arccoord{ii}(1))-varOVal];
        end
    end
    if OCMATINDIF.stateconstraint{ii}
        jumparg=freepar(OCMATINDIF.entrytimecoordinate{ii});
    end
    if OCMATINDIF.variationalcalculation
        %vswitchtimes=freepar(OCMATINDIF.vfreeswitchingtimecoordinate);
        if OCMATINDIF.stateconstraint{ii}
            vjumparg=freepar(OCMATINDIF.ventrytimecoordinate{ii});
        end
    end
    if OCMATINDIF.transversalityconditioncs{ii}
        restrans=[restrans; ...
            OCMATINDIF.bctransversalitysc(endtime,depvarb(OCMATINDIF.statecostatecoord,OCMATINDIF.arccoord{ii}(end)),jumparg(end),modelpar,OCMATINDIF.jumpid{ii}(end))];
        if OCMATINDIF.variationalcalculation
            restrans=[restrans; ...
                OCMATINDIF.variationalbctransversalitysc(endtime,depvarb(:,OCMATINDIF.arccoord{ii}),[jumparg(end) vjumparg(end)],modelpar,OCMATINDIF.jumpid{ii}(end));];
        end
    else
        restrans=[restrans; ...
            OCMATINDIF.bctransversality(endtime,depvarb(OCMATINDIF.statecostatecoord,OCMATINDIF.arccoord{ii}(end)),modelpar,arcarg(end))];
        if OCMATINDIF.variationalcalculation
            if OCMATINDIF.transversalityconditioncs{ii}
                restrans=[restrans; ...
                    OCMATINDIF.variationalbctransversality(endtime,depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,OCMATINDIF.jumpid{ii}(end));];
            else
                restrans=[restrans; ...
                    OCMATINDIF.variationalbctransversality(endtime,depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,arcarg(end));];
            end
        end
    end
    if OCMATINDIF.exogenousfunction
        resinit=[resinit; ...
            depvara(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(1))-OCMATINDIF.exogenousinitialstates(endtime,depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,arcarg)];
%         resinit=[resinit; ...
%             depvara(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(1))-OCMATINDIF.exogenousinitialstates];
    end
    if OCMATINDIF.userbc
        userbc=[userbc; ...
            OCMATINDIF.userbcfunc([OCMATINDIF.initialtime switchtimes(:).' endtime],depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,arcarg,ii)];
    end

    for arc=1:OCMATINDIF.numarc(ii)-1
        if OCMATINDIF.stateconstraint{ii} && OCMATINDIF.entryindex{ii}(arc)
            resconnec=[resconnec; ...
                OCMATINDIF.bcstateconstraint(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,jumparg(OCMATINDIF.entryindex{ii}(arc)),switchtimes,OCMATINDIF.jumpid{ii}(arc),OCMATINDIF.edge{ii},arc); ...
                OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATINDIF.edge{ii},arc)];
        else
            resconnec=[resconnec; ...
                OCMATINDIF.reset(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
                OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
        end
        if OCMATINDIF.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(arc))-depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(arc+1))];
            if OCMATINDIF.variationalcalculation
                resconnec=[resconnec; ...
                    depvarb(OCMATINDIF.variationalobjectivevaluecoord,OCMATINDIF.arccoord{ii}(arc))-depvara(OCMATINDIF.variationalobjectivevaluecoord,OCMATINDIF.arccoord{ii}(arc+1))];
            end
        end
        if OCMATINDIF.exogenousfunction
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(arc))-depvara(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(arc+1))];
        end
        if OCMATINDIF.variationalcalculation
            if OCMATINDIF.stateconstraint{ii} && OCMATINDIF.entryindex{ii}(arc)
                resconnec=[resconnec; ...
                    OCMATINDIF.variationalbcstateconstraint(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,[jumparg(OCMATINDIF.entryindex{ii}(arc)) vjumparg(OCMATINDIF.entryindex{ii}(arc))],switchtimes,OCMATINDIF.jumpid{ii}(arc),OCMATINDIF.edge{ii},arc); ...
                    OCMATINDIF.variationalguard(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATINDIF.edge{ii},arc)];
            else
                resconnec=[resconnec; ...
                    OCMATINDIF.variationalreset(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
                    OCMATINDIF.variationalguard(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
            end
        end
    end
    if OCMATINDIF.static
        O(ii)=depvarb(OCMATINDIF.static,OCMATINDIF.arccoord{ii}(end));
    else
        if OCMATINDIF.objectivevaluecalc
            O(ii)=depvarb(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(end));
        end
    end
end
for ii=1:OCMATINDIF.indifferenceorder-1
        residpt=[residpt; ...
            O(ii)-O(ii+1)];
end
res=[resinit;restarget;restrans;resconnec;residpt;userbc];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATINDIF OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
modelpar0=modelpar;
b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATINDIF.solutionindex(arc);
    modelpar=modelpar0{solutionindex};
    if OCMATINDIF.static
        modelpar(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{solutionindex});
    end
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    if OCMATINDIF.continuationtype==1
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.initialendtime+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
    end
    if ~isempty(OCMATINDIF.freetimevector)
        arctime(end)=freepar(OCMATINDIF.freetimevectorcoordinate);
    end
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
        [constr,labelS]=OCMATINDIF.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    % test jumparguments
    if 0%OCMATINDIF.stateconstraint{solutionindex}
        jumparg=freepar(OCMATINDIF.entrytimecoordinate{solutionindex});
        [constr,labelS]=OCMATINDIF.testjumpargument(arcarg,jumparg);
        violationmat=constr<-OCMATCONT.OPTIONS.admissibletol;
        if any(violationmat)
            counter=counter+1;
            cols=find(violationmat);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='jumpargument';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            %infoS(counter).constraintvalue=diffarctime(arc);
            infoS(counter).minval=min(constr);
            b=min([b infoS(counter).minval]);
        end
    end

end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
%sol=evalatmesh(tmesh,y,z);
modelpar0=modelpar;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
for arc=1:sum(OCMATINDIF.cumsumnumarc(end))
    solutionindex=OCMATINDIF.solutionindex(arc);
    modelpar=modelpar0{solutionindex};
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
    end
    if OCMATINDIF.continuationtype==1
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.initialendtime+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
    end
    diffarctime=diff(arctime);
    if solutionindex>1
        transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        transformedtimeshift=0;
    end
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(relarc)*(s(leftarcindex(arc):rightarcindex(arc))-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
end
% clear possible persistent variable
h=OCMATINDIF.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT OCBVP
idx=[];
if isempty(coeff)
    return
end
if isempty(OCMATINDIF.hitvalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
else
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

    ya = y(:,OCBVP.Lidx);
    yb = y(:,OCBVP.Ridx);

    if OCMATINDIF.static
        for ii=1:OCMATINDIF.indifferenceorder
            modelpar(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{ii});
            modelpartmp{ii}=modelpar;
        end
    else
        modelpartmp=modelpar;
    end
    fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
    fprintf(1,' Difference to hitvalue: %g\n',OCMATINDIF.hitvaluefunc(t,ya,yb,modelpartmp,OCMATCONT.HE.arcarg,OCMATINDIF.arccoord,OCMATINDIF.hitvalue));
end

%fprintf(1,' Continuation parameter: %g\n',coeff(end));
%fprintf(1,' Continuation parameter: %g\n',coeff(end));

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
global OCMATCONT OCMATINDIF OCBVP

failed=[];
if isempty(OCMATINDIF.hitvalue)
    out=1-coeff(end);
else
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

    ya = y(:,OCBVP.Lidx);
    yb = y(:,OCBVP.Ridx);

    if OCMATINDIF.static
        for ii=1:OCMATINDIF.indifferenceorder
            modelpar(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{ii});
            modelpartmp{ii}=modelpar;
        end
    else
        modelpartmp=modelpar{1};
    end
    out=OCMATINDIF.hitvaluefunc(t,ya,yb,modelpartmp,OCMATCONT.HE.arcarg,OCMATINDIF.arccoord,OCMATINDIF.hitvalue);
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,modelpar);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[];
out.timehorizon=[];
for ii=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.continuationtype==1
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' OCMATINDIF.initialendtime+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' OCMATINDIF.endtime];
    end
    if ~isempty(OCMATINDIF.freetimevector)
        arctime(end)=freepar(OCMATINDIF.freetimevectorcoordinate);
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
    if OCMATINDIF.static
        out.modelparameter{ii}(OCMATINDIF.staticparameterindex)=freepar(OCMATINDIF.staticparametercoordinate{ii});
        out.modelparameter{ii}=modelpar{ii};
    else
        out.modelparameter{ii}=modelpar{ii};
    end

end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATINDIF.initialtime;
%out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution4ft';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.initcoord=OCMATINDIF.initcoord;
out.solverinfo.solutionindex=OCMATINDIF.solutionindex;
out.solverinfo.initialstateindex=OCMATINDIF.initialstateindex;
out.solverinfo.indifferenceorder=OCMATINDIF.indifferenceorder;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.statecostatecoord=OCMATINDIF.statecostatecoord;
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.objectivevaluecoord=OCMATINDIF.objectivevaluecoord;
out.solverinfo.parameterindex=OCMATINDIF.parameterindex;
out.solverinfo.freeparameterindex=OCMATINDIF.freeparameterindex;
out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.solverinfo.continuationtype=OCMATINDIF.continuationtype;
out.solverinfo.continuationindex=OCMATINDIF.continuationindex;
out.solverinfo.continuationvector=OCMATINDIF.continuationvector;
out.solverinfo.continuationcoordinate=OCMATINDIF.continuationcoordinate;
out.solverinfo.stateconstraint=OCMATINDIF.stateconstraint;
out.solverinfo.static=OCMATINDIF.static;
if any([OCMATINDIF.stateconstraint{:}])
    out.solverinfo.entrytimecoordinate=OCMATINDIF.entrytimecoordinate;
    out.solverinfo.transversalityconditioncs=OCMATINDIF.transversalityconditioncs;
    out.solverinfo.entryindex=OCMATINDIF.entryindex;
    out.solverinfo.jumpid=OCMATINDIF.jumpid;
    out.solverinfo.jumpargument=OCMATINDIF.jumpargument;
    if OCMATINDIF.variationalcalculation
        %out.solverinfo.vfreeparametercoordinate=OCMATINDIF.vfreeparametercoordinate;
        for ii=1:OCMATINDIF.indifferenceorder
            out.solverinfo.vjumpargument{ii}=freepar(OCMATINDIF.ventrytimecoordinate{ii});
        end
    end
end
if OCMATINDIF.variationalcalculation
    out.solverinfo.variationalobjectivevaluecoord=OCMATINDIF.variationalobjectivevaluecoord;
    for ii=1:OCMATINDIF.indifferenceorder
        out.solverinfo.vfreetime{ii}=freepar(OCMATINDIF.vfreetimecoord{ii});
    end
end
if OCMATINDIF.static
    out.solverinfo.staticparameterindex=OCMATINDIF.staticparameterindex;
    out.solverinfo.staticparametercoordinate=OCMATINDIF.staticparametercoordinate;
end

switch OCMATINDIF.continuationtype
    case 0
        out.solverinfo.initialstate=OCMATINDIF.initialstate;
    case 1
        out.solverinfo.initialendtime=OCMATINDIF.initialendtime;
    case 2
        out.solverinfo.initialparameter=OCMATINDIF.initialparameter;
end
out.solverinfo.numarc=OCMATINDIF.numarc;
out.solverinfo.arccoord=OCMATINDIF.arccoord;
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

modelpar=OCMATINDIF.modelparameter;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
for ii=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.continuationtype==2
        modelpar{ii}(OCMATINDIF.continuationindex)=OCMATINDIF.initialparameter+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
    end
    if OCMATINDIF.freeparametervector
        if isempty(OCMATINDIF.freeinitialparametervalue)
            modelpar{ii}(OCMATINDIF.freeparameterindex)=freepar(OCMATINDIF.freeparametervectorcoordinate);
        else
            modelpar{ii}(OCMATINDIF.freeparameterindex)=OCMATINDIF.freeinitialparametervalue+freepar(OCMATINDIF.freeparametervectorcoordinate)*OCMATINDIF.freeparameterdirection;
        end
    end
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
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4ft'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4ft'],'sout','bvpout')
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
        mbcidx=find(~diff(tmesh));
        Lidx = [1, mbcidx+1];
        Ridx = [mbcidx, length(tmesh)];
        mbcidx=find(~diff(tmeshnew));
        Lidxnew = [1, mbcidx+1];
        Ridxnew = [mbcidx, length(tmeshnew)];
        warning('off','MATLAB:deval:NonuniqueSolution');
        ynew=zeros(size(sol.y,1),length(tmeshnew));
        solpart.solver=sol.solver;
        for ii=1:length(Lidx)
            solpart.y=sol.y(:,Lidx(ii):Ridx(ii));
            solpart.yp=sol.yp(:,Lidx(ii):Ridx(ii));
            solpart.x=sol.x(Lidx(ii):Ridx(ii));
            ynew(:,Lidxnew(ii):Ridxnew(ii))=devalbvpoc(solpart,tmeshnew(Lidxnew(ii):Ridxnew(ii)));
        end
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