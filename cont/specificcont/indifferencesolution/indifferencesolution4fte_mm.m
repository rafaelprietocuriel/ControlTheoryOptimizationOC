function out=indifferencesolution4fte_mm()

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

if OCMATINDIF.continuationtype==2
    for ii=1:OCMATINDIF.totalnumberofparts
        if ~isempty(OCMATINDIF.initialparameter{ii})
            modelpar{ii}(OCMATINDIF.continuationindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
        end
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for ii=1:OCMATINDIF.totalnumberofparts
         modelpar{ii}(OCMATINDIF.freeparameterindex{ii})=freepar(OCMATINDIF.freeparametervectorcoordinate{ii});
    end
end

funcindex=OCMATINDIF.funcindex(arc);
indifferenceindex=OCMATINDIF.indifferenceindex(arc);
relarcindex=OCMATINDIF.relarcindex(arc);
transformedtimeshift=OCMATINDIF.timeshift(indifferenceindex);
arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{indifferenceindex});
arctime(OCMATINDIF.timepointsfreeindex{indifferenceindex})=freepar(OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex});
arctime(OCMATINDIF.timepointsoptimizedindex{indifferenceindex})=freepar(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex});
if OCMATINDIF.continuationtype==1
    %arctime(OCMATINDIF.continuationindex(1))=OCMATINDIF.initialarcinterval(OCMATINDIF.continuationindex(1))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(1);
    arctime(OCMATINDIF.continuationindex(indifferenceindex))=arctime(OCMATINDIF.continuationindex(indifferenceindex))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(indifferenceindex);
end
if ~isempty(OCMATINDIF.freetimevector)
    arctime(OCMATINDIF.freetimevector(indifferenceindex))=freepar(OCMATINDIF.freetimevectorindex);
end

ct=arctime(OCMATINDIF.connectiontimeindex{indifferenceindex});

diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarcindex)*(s-transformedtimeshift)+(arctime(relarcindex)-diffarctime(relarcindex)*(relarcindex-1));

dtds=diffarctime(relarcindex);
dxdt=dtds*OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
if OCMATINDIF.explicitconnectiontime
    dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=dtds*OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
end

if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
end
if OCMATINDIF.exogenousfunction
    dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=dtds*OCMATINDIF.exogenousdynamics{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
partindex=OCMATINDIF.relpartindex(arc);
funcindex=OCMATINDIF.funcindex(arc);
indifferenceindex=OCMATINDIF.indifferenceindex(arc);
relarcindex=OCMATINDIF.relarcindex(arc);
transformedtimeshift=OCMATINDIF.timeshift(indifferenceindex);
arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{indifferenceindex});
arctime(OCMATINDIF.timepointsfreeindex{indifferenceindex})=freepar(OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex});
arctime(OCMATINDIF.timepointsoptimizedindex{indifferenceindex})=freepar(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex});
if OCMATINDIF.continuationtype==1
    arctime(OCMATINDIF.continuationindex(indifferenceindex))=arctime(OCMATINDIF.continuationindex(indifferenceindex))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(indifferenceindex);
end
if ~isempty(OCMATINDIF.freetimevector)
    arctime(OCMATINDIF.freetimevector(indifferenceindex))=freepar(OCMATINDIF.freetimevectorindex);
end

ct=arctime(OCMATINDIF.connectiontimeindex{indifferenceindex});

diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(relarcindex)*(s-transformedtimeshift)+(arctime(relarcindex)-diffarctime(relarcindex)*(relarcindex-1));
dtds=diffarctime(relarcindex(ones(1,numel(s))));
J=OCMATINDIF.canonicalsystemjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
J=dtds*J;
if OCMATINDIF.explicitconnectiontime
    J=[J;dtds*OCMATINDIF.dhamiltoniandctjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct)];
end
if OCMATINDIF.objectivevaluecalc
    J=[J;dtds*OCMATINDIF.objectivefunctionjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct)];
end
if OCMATINDIF.exogenousfunction
    J=[J;dtds*OCMATINDIF.exogenousjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct)];
end
J=[J,OCMATINDIF.Jext];
% if OCMATINDIF.explicitconnectiontime
%     J=[J; ...
%         dtds*OCMATINDIF.dhamiltoniandctjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct)];
%     J=[J OCMATINDIF.zeros2nxct];
% end
% if OCMATINDIF.objectivevaluecalc
%         J=[J; ...
%             [dtds*OCMATINDIF.objectivefunctionjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct) OCMATINDIF.zerosct]];
%         J=[J OCMATINDIF.zeros2npo];
% end
% if OCMATINDIF.exogenousfunction
%         J=[J; ...
%             [dtds*OCMATINDIF.exogenousjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct) OCMATINDIF.zerosctpo]];
%         J=[J OCMATINDIF.zeros2npef];
% end
%Jpar=OCMATINDIF.zerosnpctxpar;
Jpar=OCMATINDIF.Jpar;
if OCMATINDIF.timepointscharacterization{indifferenceindex}(relarcindex+1)==1
    dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    if OCMATINDIF.explicitconnectiontime
        dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    Jpar(:,OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex}(relarcindex+1==OCMATINDIF.timepointsfreeindex{indifferenceindex}))=dxdt;
end
if OCMATINDIF.timepointscharacterization{indifferenceindex}(relarcindex)==1
    dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    if OCMATINDIF.explicitconnectiontime
        dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    Jpar(:,OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex}(relarcindex==OCMATINDIF.timepointsfreeindex{indifferenceindex}))=-dxdt;
end
if OCMATINDIF.timepointscharacterization{indifferenceindex}(relarcindex+1)==-1
    dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    if OCMATINDIF.explicitconnectiontime
        dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    if OCMATINDIF.exogenousfunction
        dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
    end
    Jpar(:,OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex}(relarcindex+1==OCMATINDIF.timepointsoptimizedindex{indifferenceindex}))=dxdt+Jt;
end
if OCMATINDIF.timepointscharacterization{indifferenceindex}(relarcindex)==-1
    dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunctionderivativetime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    if OCMATINDIF.explicitconnectiontime
        dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        tmp=dtds*OCMATINDIF.derivativeconnectiontime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jt(OCMATINDIF.statecostatecoord,:)=tmp(OCMATINDIF.statecostatecoord,relarcindex);
        tmp=dtds*OCMATINDIF.d2hamiltoniandct2{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jt(OCMATINDIF.dhamiltoniandctcoord,:)=tmp(relarcindex);
        if OCMATINDIF.objectivevaluecalc
            tmp=dtds*OCMATINDIF.objectivefunctionderivativeconnectiontime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
            Jt(OCMATINDIF.objectivevaluecoord,:)=Jt(OCMATINDIF.objectivevaluecoord,:)+tmp(relarcindex);
        end
    end
    if OCMATINDIF.exogenousfunction
        dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
    end
    Jpar(:,OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex}(relarcindex==OCMATINDIF.timepointsoptimizedindex{indifferenceindex}))=-dxdt+Jt;
end


%%%%%
if OCMATINDIF.continuationtype==1
    if  OCMATINDIF.timepointscharacterization{indifferenceindex}(relarcindex)==2
        dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        if OCMATINDIF.objectivevaluecalc
            dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
            Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        end
        if OCMATINDIF.explicitconnectiontime
            dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        end
        if OCMATINDIF.exogenousfunction
            dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
            Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
        end
        idx=find(OCMATINDIF.continuationindex==relarcindex);
        Jpar(:,OCMATINDIF.continuationcoordinate)=Jpar(:,OCMATINDIF.continuationcoordinate)-OCMATINDIF.continuationvector(idx)*dxdt;
    end
    if  OCMATINDIF.timepointscharacterization{indifferenceindex}(relarcindex+1)==2
        dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        if OCMATINDIF.objectivevaluecalc
            dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
            Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        end
        if OCMATINDIF.explicitconnectiontime
            dxdt(OCMATINDIF.dhamiltoniandctcoord,:)=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        end
        if OCMATINDIF.exogenousfunction
            dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
            Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
        end
        idx=find(OCMATINDIF.continuationindex==relarcindex+1);
        Jpar(:,OCMATINDIF.continuationcoordinate)=Jpar(:,OCMATINDIF.continuationcoordinate)+OCMATINDIF.continuationvector(idx)*dxdt;
    end
    %     if arc==OCMATCONT.HE.numarc  && any(OCMATINDIF.continuationindex==OCMATCONT.HE.numarc+1)
    %         dxdt=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    %         if OCMATINDIF.objectivevaluecalc
    %             dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    %             Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    %         end
    %         idx=find(OCMATINDIF.continuationindex==OCMATCONT.HE.numarc+1);
    %         Jpar(:,OCMATINDIF.continuationcoordinate)=Jpar(:,OCMATINDIF.continuationcoordinate)-OCMATINDIF.continuationvector(idx)*dxdt;
    %     end
end
if ~isempty(OCMATINDIF.freeparametervector)
    Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    for jj=1:OCMATINDIF.numberoffreeparameter
        if ~isempty(OCMATINDIF.parameterindex{jj,funcindex})
            Jpar(OCMATINDIF.statecostatecoord,OCMATINDIF.freeparametervectorcoordinate(jj))=Jmodelpar(:,OCMATINDIF.parameterindex{jj,funcindex});
        end
    end
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        for jj=1:OCMATINDIF.numberoffreeparameter
            if ~isempty(OCMATINDIF.parameterindex{jj,funcindex})
                Jpar(OCMATINDIF.objectivevaluecoord,OCMATINDIF.freeparametervectorcoordinate(jj))=Jobjmodelpar(:,OCMATINDIF.parameterindex{jj,funcindex});
            end
        end
    end
    if OCMATINDIF.explicitconnectiontime
        Jdhdctmodelpar=dtds*OCMATINDIF.dhamiltoniandctparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        for jj=1:OCMATINDIF.numberoffreeparameter
            if ~isempty(OCMATINDIF.parameterindex{jj,funcindex})
                Jpar(OCMATINDIF.dhamiltoniandctcoord,OCMATINDIF.freeparametervectorcoordinate(jj))=Jdhdctmodelpar(:,OCMATINDIF.parameterindex{jj,funcindex});
            end
        end
    end
    if OCMATINDIF.exogenousfunction
        Jexfmodelpar=dtds*OCMATINDIF.exogenousparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        for jj=1:OCMATINDIF.numberoffreeparameter
            if ~isempty(OCMATINDIF.parameterindex{jj,funcindex})
                Jpar(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.freeparametervectorcoordinate(jj))=Jexfmodelpar(:,OCMATINDIF.parameterindex{jj,funcindex});
            end
        end
    end
end
if ~isempty(OCMATINDIF.freetimevector) && OCMATINDIF.freetimevector(indifferenceindex)==relarcindex+1
    Jmodelpar=OCMATINDIF.canonicalsystem{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    Jpar(OCMATINDIF.statecostatecoord,OCMATINDIF.freetimevectorindex)=Jmodelpar(OCMATINDIF.statecostatecoord);
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=OCMATINDIF.objectivefunction{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jpar(OCMATINDIF.objectivevaluecoord,OCMATINDIF.freetimevectorindex)=Jobjmodelpar;
    end
    if OCMATINDIF.explicitconnectiontime
        Jdhdctmodelpar=OCMATINDIF.dhamiltoniandct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jpar(OCMATINDIF.dhamiltoniandctcoord,OCMATINDIF.freetimevectorindex)=Jdhdctmodelpar;
    end
    if OCMATINDIF.exogenousfunction
        Jexfmodelpar=OCMATINDIF.exogenousdynamics{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        Jpar(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.freetimevectorindex)=Jexfmodelpar;
    end
end

if OCMATINDIF.continuationtype==2
    Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);

    if length(OCMATINDIF.continuationindex)==1
        Jpar(OCMATINDIF.statecostatecoord,end)=OCMATINDIF.continuationvector{funcindex}*Jmodelpar(:,OCMATINDIF.continuationindex);
    else
        Jpar(OCMATINDIF.statecostatecoord,end)=sum(repmat(OCMATINDIF.continuationvector{funcindex},size(Jmodelpar(:,OCMATINDIF.continuationindex{funcindex}),1),1).*Jmodelpar(:,OCMATINDIF.continuationindex{funcindex}),2);
    end
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        if length(OCMATINDIF.continuationindex)==1
            Jpar(OCMATINDIF.objectivevaluecoord,end)=OCMATINDIF.continuationvector{funcindex}*Jobjmodelpar(:,OCMATINDIF.continuationindex);
        else
            Jpar(OCMATINDIF.objectivevaluecoord,end)=sum(repmat(OCMATINDIF.continuationvector{funcindex},size(Jobjmodelpar(:,OCMATINDIF.continuationindex{funcindex}),1),1).*Jobjmodelpar(:,OCMATINDIF.continuationindex{funcindex}),2);
        end
    end
    if OCMATINDIF.explicitconnectiontime
        Jdhdctmodelpar=dtds*OCMATINDIF.dhamiltoniandctparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        if length(OCMATINDIF.continuationindex)==1
            Jpar(OCMATINDIF.dhamiltoniandctcoord,end)=OCMATINDIF.continuationvector{funcindex}*Jdhdctmodelpar(:,OCMATINDIF.continuationindex);
        else
            Jpar(OCMATINDIF.dhamiltoniandctcoord,end)=sum(repmat(OCMATINDIF.continuationvector{funcindex},size(Jdhdctmodelpar(:,OCMATINDIF.continuationindex{funcindex}),1),1).*Jdhdctmodelpar(:,OCMATINDIF.continuationindex{funcindex}),2);
        end
    end
    if OCMATINDIF.exogenousfunction
        Jexfmodelpar=dtds*OCMATINDIF.exogenousparameterjacobian{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
        if length(OCMATINDIF.continuationindex)==1
            Jpar(OCMATINDIF.exogenousdynamicscoordinate,end)=OCMATINDIF.continuationvector{funcindex}*Jexfmodelpar(:,OCMATINDIF.continuationindex);
        else
            Jpar(OCMATINDIF.exogenousdynamicscoordinate,end)=sum(repmat(OCMATINDIF.continuationvector{funcindex},size(Jexfmodelpar(:,OCMATINDIF.continuationindex{funcindex}),1),1).*Jobjmodelpar(:,OCMATINDIF.continuationindex{funcindex}),2);
        end
    end
end
if OCMATINDIF.explicitconnectiontime && ~isempty(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex})
    d2xdtdc=dtds*OCMATINDIF.derivativeconnectiontime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    if OCMATINDIF.objectivevaluecalc
        d2xdtdc(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunctionderivativeconnectiontime{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    if OCMATINDIF.exogenousfunction
        d2xdtdc(OCMATINDIF.exogenousdynamicscoordinate,:)=dtds*OCMATINDIF.exogenousdct{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    end
    d2xdtdc(OCMATINDIF.dhamiltoniandctcoord,:)=dtds*OCMATINDIF.d2hamiltoniandct2{funcindex}(t,depvar,modelpar{funcindex},arcarg,ct);
    Jpar(:,OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex})=Jpar(:,OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex})+d2xdtdc(:,OCMATINDIF.timepointsoptimizedindex{indifferenceindex}-1);
end


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT
resconnec=[];
residpt=[];
resinit=[];
restrans=[];
restarget=[];
respartconnect=[];
if OCMATINDIF.continuationtype==2
    for ii=1:OCMATINDIF.totalnumberofparts
        if ~isempty(OCMATINDIF.initialparameter{ii})
            modelpar{ii}(OCMATINDIF.continuationindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
        end
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for ii=1:OCMATINDIF.totalnumberofparts
         modelpar{ii}(OCMATINDIF.freeparameterindex{ii})=freepar(OCMATINDIF.freeparametervectorcoordinate{ii});
    end
end


O=zeros(1,OCMATINDIF.indifferenceorder);
if OCMATINDIF.explicitconnectiontime
    dHdctval=zeros(OCMATINDIF.connectiontimenumber,OCMATINDIF.indifferenceorder);
end
initialstate=OCMATINDIF.initialstate(:,1);
if OCMATINDIF.continuationtype==0
    initialstate=initialstate+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
end
if ~isempty(OCMATINDIF.freestatevector)
    for jj=1:length(OCMATINDIF.freestatevectorcoordinate)
        initialstate=initialstate+freepar(OCMATINDIF.freestatevectorcoordinate(jj))*OCMATINDIF.freestatevector(:,jj);
    end
end
for arc=1:OCMATINDIF.totalnumberofarcs
    part=OCMATINDIF.relpartindex(arc);
    relarcindex=OCMATINDIF.relarcindex(arc);
    funcindex=OCMATINDIF.funcindex(arc);
    indifferenceindex=OCMATINDIF.indifferenceindex(arc);
    %depvarpartcoord=OCMATINDIF.indiffence2partcoordinate{arc};
    arcarg=OCMATCONT.HE.arcarg(arc);
    if relarcindex==1
        OVal=0;
        depvarcoor=OCMATINDIF.relcoord{indifferenceindex};
        arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{indifferenceindex});
        arctime(OCMATINDIF.timepointsfreeindex{indifferenceindex})=freepar(OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex});
        arctime(OCMATINDIF.timepointsoptimizedindex{indifferenceindex})=freepar(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex});
        if OCMATINDIF.continuationtype==1
            arctime(OCMATINDIF.continuationindex(indifferenceindex))=arctime(OCMATINDIF.continuationindex(indifferenceindex))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(indifferenceindex);
        end
        if ~isempty(OCMATINDIF.freetimevector)
            arctime(OCMATINDIF.freetimevector(indifferenceindex))=freepar(OCMATINDIF.freetimevectorindex);
        end
        ct=arctime(OCMATINDIF.connectiontimeindex{indifferenceindex});
        if ~isempty(OCMATINDIF.userbc) && OCMATINDIF.userbc{indifferenceindex}
            restarget=[restarget; ...
                OCMATINDIF.userfunctionbc{arc}(depvara,depvarb,modelpar,ct,OCMATCONT.HE.arcarg,arc,part)];
        end

        switchtimes=arctime(2:end-1);
        relarcargument=OCMATINDIF.relarcargument{indifferenceindex};
        if indifferenceindex==1%<OCMATINDIF.indifferenceorder
            if ~isempty(OCMATINDIF.fixinitialstate)
                resinit=[resinit; ...
                    depvara(OCMATINDIF.fixinitialstatecoordinate,1)-OCMATINDIF.initialstate(OCMATINDIF.fixinitialstatecoordinate,1)];
            end
            if OCMATINDIF.continuationtype==0
                resinit=[resinit; ...
                    depvara(OCMATINDIF.continuationindex,1)-initialstate];
            end
        end
    end
    if OCMATINDIF.changepart(arc)==-1 % change solution path, last arc of last solution paths
        if arc<OCMATINDIF.totalnumberofarcs
            depvarcoord1=OCMATINDIF.relcoord{indifferenceindex+1};
            resinit=[resinit; ...
                depvara(OCMATINDIF.statecoordinate,depvarcoor(:,1))-depvara(OCMATINDIF.statecoordinate,depvarcoord1(:,1))];
        end
        restrans=[restrans; ...
            OCMATINDIF.bctransversality{funcindex}(arctime(end),depvarb(OCMATINDIF.statecostatecoord,depvarcoor(:,end)),modelpar{funcindex},arcarg(end),ct)];
        if OCMATINDIF.objectivevaluecalc
            OVal=OVal+OCMATINDIF.salvagevalue{funcindex}(arctime(end),depvarb(OCMATINDIF.statecostatecoord,depvarcoor(:,end)),modelpar{funcindex},arcarg(end),ct);
            resinit=[resinit; ...
                depvara(OCMATINDIF.objectivevaluecoord,depvarcoor(1))-OVal];
            O(indifferenceindex)=depvarb(OCMATINDIF.objectivevaluecoord,depvarcoor(:,end));
        end
        if OCMATINDIF.explicitconnectiontime
            dHdctval(:,indifferenceindex)=dHdctval(:,indifferenceindex)+OCMATINDIF.dsalvagedct{funcindex}(arctime(relarcindex+1),depvarb(OCMATINDIF.statecostatecoord,relarcindex),modelpar{funcindex},OCMATCONT.HE.arcarg(arc),ct);
            resinit=[resinit; ...
                depvara(OCMATINDIF.dhamiltoniandctcoord,depvarcoor(1))-dHdctval(:,indifferenceindex)];
        end
        if OCMATINDIF.exogenousfunction
            resinit=[resinit; ...
                depvara(OCMATINDIF.exogenousdynamicscoordinate,depvarcoor(1))-OCMATINDIF.exogenousinitialstates];
        end

    elseif OCMATINDIF.changepart(arc)==1 % change stage within a solution path
        respartconnect=[respartconnect; ...
            OCMATINDIF.bcconnectingparts{funcindex}(depvara(:,depvarcoor),depvarb(:,depvarcoor),modelpar(OCMATINDIF.modelparameterindex{indifferenceindex}),ct,relarcargument,relarcindex)];
        if any(relarcindex==OCMATINDIF.timepointsoptimizedindex{indifferenceindex}-1)
            respartconnect=[respartconnect; ...
                OCMATINDIF.bcoptimalconnectingparts{funcindex}(depvara(:,depvarcoor),depvarb(:,depvarcoor),modelpar(OCMATINDIF.modelparameterindex{indifferenceindex}),ct,relarcargument,relarcindex,part)];
        end
        if OCMATINDIF.objectivevaluecalc
            %OVal=OVal+OCMATINDIF.salvagevalue{funcindex}(arctime(arc),depvarb(OCMATINDIF.statecostatecoord,relarcindex),modelpar{funcindex},OCMATCONT.HE.arcarg(arc),ct);
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.objectivevaluecoord,depvarcoor(relarcindex))-depvara(OCMATINDIF.objectivevaluecoord,depvarcoor(relarcindex+1))];
        end
        if OCMATINDIF.exogenousfunction
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.exogenousdynamicscoordinate,depvarcoor(relarcindex))-depvara(OCMATINDIF.exogenousdynamicscoordinate,depvarcoor(relarcindex+1))];
        end
        if OCMATINDIF.explicitconnectiontime
            dHdctval(:,indifferenceindex)=dHdctval(:,indifferenceindex)+OCMATINDIF.dsalvagedct{funcindex}(arctime(relarcindex+1),depvarb(OCMATINDIF.statecostatecoord,relarcindex),modelpar{funcindex},OCMATCONT.HE.arcarg(arc),ct);
            respartconnect=[respartconnect; ...
                depvara(OCMATINDIF.dhamiltoniandctcoord,depvarcoor(relarcindex+1))-depvarb(OCMATINDIF.dhamiltoniandctcoord,depvarcoor(relarcindex))];
        end
    elseif ~OCMATINDIF.changepart(arc) % change arc within a part of a solution path
        resconnec=[resconnec; ...
            OCMATINDIF.reset{arc}(depvara(:,depvarcoor),depvarb(:,depvarcoor),modelpar{arc},switchtimes{indifferenceindex}(OCMATINDIF.parts2timepointsfreeindex{arc}),arcarg,OCMATINDIF.edge{arc},arc); ...
            OCMATINDIF.guard{arc}(depvara(:,depvarcoor),depvarb(:,depvarcoor),modelpar{arc},switchtimes{indifferenceindex},arcarg,OCMATINDIF.edge{arc},arc)];
        if OCMATINDIF.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.objectivevaluecoord,depvarcoor(relarcindex))-depvara(OCMATINDIF.objectivevaluecoord,depvarcoor(relarcindex+1))];
        end
        if OCMATINDIF.exogenousfunction
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.exogenousdynamicscoordinate,depvarcoor(relarcindex))-depvara(OCMATINDIF.exogenousdynamicscoordinate,depvarcoor(relarcindex+1))];
        end
        if OCMATINDIF.explicitconnectiontime
            respartconnect=[respartconnect; ...
                depvara(OCMATINDIF.dhamiltoniandctcoord,arc+1)-depvarb(OCMATINDIF.dhamiltoniandctcoord,arc)];
        end
    end
end
for arc=1:OCMATINDIF.indifferenceorder-1
    residpt=[residpt; ...
        O(arc)-O(arc+1)];
end

res=[resinit;restarget;restrans;resconnec;residpt;respartconnect];
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

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    partindex=OCMATINDIF.relpartindex(arc);
    funcindex=OCMATINDIF.funcindex(arc);
    indifferenceindex=OCMATINDIF.indifferenceindex(arc);
    relarcindex=OCMATINDIF.relarcindex(arc);
    transformedtimeshift=OCMATINDIF.timeshift(indifferenceindex);
    arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{indifferenceindex});
    arctime(OCMATINDIF.timepointsfreeindex{indifferenceindex})=freepar(OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex});
    arctime(OCMATINDIF.timepointsoptimizedindex{indifferenceindex})=freepar(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex});
    if OCMATINDIF.continuationtype==1
        arctime(OCMATINDIF.continuationindex(indifferenceindex))=arctime(OCMATINDIF.continuationindex(indifferenceindex))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(indifferenceindex);
    end
    if ~isempty(OCMATINDIF.freetimevector)
        arctime(OCMATINDIF.freetimevector(indifferenceindex))=freepar(OCMATINDIF.freetimevectorindex);
    end
    ct=arctime(OCMATINDIF.connectiontimeindex{indifferenceindex});
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarcindex)*(s-transformedtimeshift)+(arctime(relarcindex)-diffarctime(relarcindex)*(relarcindex-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATINDIF.testadmissibility{funcindex}(t,sol.y(idx,:),modelpar{indifferenceindex},arcarg,ct);
    else
        eqcoord=domainddata{indifferenceindex}(arcindex).eqcoord;
        [constr labelS]=OCMATINDIF.testadmissibility{funcindex}(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar{indifferenceindex},arcarg,ct);
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
    violationmat=diffarctime(relarcindex)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
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
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
for indifferenceindex=1:OCMATINDIF.indifferenceorder
    arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{indifferenceindex});
    arctime(OCMATINDIF.timepointsfreeindex{indifferenceindex})=freepar(OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex});
    arctime(OCMATINDIF.timepointsoptimizedindex{indifferenceindex})=freepar(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex});
    if OCMATINDIF.continuationtype==1
        arctime(OCMATINDIF.continuationindex(indifferenceindex))=arctime(OCMATINDIF.continuationindex(indifferenceindex))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(indifferenceindex);
    end
    if ~isempty(OCMATINDIF.freetimevector)
        arctime(OCMATINDIF.freetimevector(indifferenceindex))=freepar(OCMATINDIF.freetimevectorindex);
    end
    diffarctime=diff(arctime);
    ctr=0;
    for arc=OCMATINDIF.relcoord{indifferenceindex}
        ctr=ctr+1;
        t(leftarcindex(arc):rightarcindex(arc))=diffarctime(ctr)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(ctr)-diffarctime(ctr)*(arc-1));
    end

end
h=OCMATINDIF.plotcontinuation{1}(t,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
if ~isempty(OCMATINDIF.freetargetvalue)
    [t,y,z,freepar]=drearr(tmesh,coeff,tangent);
    fprintf(1,' Difference to target value: %g\n',OCMATINDIF.freetargetvalue-freepar(end-1));
else
    fprintf(1,' Continuation parameter: %g\n',coeff(end));
end



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
global OCMATCONT OCMATINDIF

failed=[];
out=[];
for ii=id
    switch ii
        case 1
            if ~isempty(OCMATINDIF.freetargetvalue)
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                out=OCMATINDIF.freetargetvalue-freepar(end-1);
            else
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(end);
            end
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
out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
    OCMATCONT.HE.TIMEDDATA.rightarcindex];
out.arcinterval=[];
out.timehorizon=[];
for indifferenceindex=1:OCMATINDIF.indifferenceorder
    arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{indifferenceindex});
    arctime(OCMATINDIF.timepointsfreeindex{indifferenceindex})=freepar(OCMATINDIF.timepointsfree4parcoordinate{indifferenceindex});
    arctime(OCMATINDIF.timepointsoptimizedindex{indifferenceindex})=freepar(OCMATINDIF.timepointsoptimize4parcoordinate{indifferenceindex});
    if OCMATINDIF.continuationtype==1
        arctime(OCMATINDIF.continuationindex(indifferenceindex))=arctime(OCMATINDIF.continuationindex(indifferenceindex))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(indifferenceindex);
    end
    if ~isempty(OCMATINDIF.freetimevector)
        arctime(OCMATINDIF.freetimevector(indifferenceindex))=freepar(OCMATINDIF.freetimevectorindex);
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
end
if OCMATINDIF.continuationtype==2
    for ii=1:OCMATINDIF.totalnumberofparts
        if ~isempty(OCMATINDIF.initialparameter{ii})
            modelpar{ii}(OCMATINDIF.continuationindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
        end
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for ii=1:OCMATINDIF.totalnumberofparts
         modelpar{ii}(OCMATINDIF.freeparameterindex{ii})=freepar(OCMATINDIF.freeparametervectorcoordinate{ii});
    end
end


out.arcarg=OCMATCONT.HE.arcarg;
out.x0=arctime(1);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution4fte_mm';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.indifferenceindex=OCMATINDIF.indifferenceindex;
out.solverinfo.indifferenceorder=OCMATINDIF.indifferenceorder;
out.solverinfo.timepointsfree4parcoordinate=OCMATINDIF.timepointsfree4parcoordinate;
out.solverinfo.statecostatecoord=OCMATINDIF.statecostatecoord;
out.solverinfo.partindex=OCMATINDIF.relpartindex;
out.solverinfo.funcindex=OCMATINDIF.funcindex;
out.solverinfo.relarcindex=OCMATINDIF.relarcindex;
out.solverinfo.transformedtimeshift=OCMATINDIF.timeshift;
out.solverinfo.freetimevector=OCMATINDIF.freetimevector;
out.solverinfo.objectivevaluecoord=OCMATINDIF.objectivevaluecoord;
out.solverinfo.dhamiltoniandctcoord=OCMATINDIF.dhamiltoniandctcoord;
out.solverinfo.exogenousdynamicscoordinate=OCMATINDIF.exogenousdynamicscoordinate;
out.solverinfo.userbc=OCMATINDIF.userbc;
out.solverinfo.continuationtype=OCMATINDIF.continuationtype;

if OCMATINDIF.continuationtype==2
    out.solverinfo.continuationindex=OCMATINDIF.continuationindex;
end
out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;

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

modelpar=OCMATINDIF.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
if OCMATINDIF.continuationtype==2
    for ii=1:OCMATINDIF.totalnumberofparts
        if ~isempty(OCMATINDIF.initialparameter{ii})
            modelpar{ii}(OCMATINDIF.continuationindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
        end
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for ii=1:OCMATINDIF.totalnumberofparts
        modelpar{ii}(OCMATINDIF.freeparameterindex{ii})=freepar(OCMATINDIF.freeparametervectorcoordinate{ii});
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
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4fte_mm'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4fte_mm'],'sout','bvpout')
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