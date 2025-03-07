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

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
arctime=OCMATFTE.initialtimeinterval;
arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
end
solutionindex=OCMATFTE.solutionindex(arc);
if OCMATFTE.freeparameter
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=freepar(OCMATFTE.parametercoord);
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=OCMATFTE.initialparameter(ii,:)+freepar(end)*OCMATFTE.continuationvector(ii,:);
    end
end

diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));

dtds=diffarctime(arc);
dxdt=dtds*OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
if OCMATFTE.objectivevaluecalc
    dxdt(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
arctime=OCMATFTE.initialtimeinterval;
arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
end
solutionindex=OCMATFTE.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATFTE.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc(ones(1,numel(s))));
J=OCMATFTE.canonicalsystemjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
J=dtds*J;
if OCMATFTE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATFTE.objectivefunctionjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).numeq,OCMATCONT.HE.numparameter);
if any(arc==OCMATFTE.switchtimeindex)
    dxdt=OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
    Jpar(:,OCMATFTE.switchtimecoordinate(arc==OCMATFTE.switchtimeindex))=-dxdt;
end
if any(arc==OCMATFTE.switchtimeindex-1)
    dxdt=OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
    Jpar(:,OCMATFTE.switchtimecoordinate(arc==OCMATFTE.switchtimeindex-1))=dxdt;
end
if any(arc==OCMATFTE.optimalparttimeindex)
    idx=find(arc==OCMATFTE.optimalparttimeindex);
    dxdt=OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
    Jpar(:,OCMATFTE.optimalpartcoordinate(idx))=-dxdt+Jt;
end
if any(arc==OCMATFTE.optimalparttimeindex-1)
    idx=find(arc==OCMATFTE.optimalparttimeindex-1);
    dxdt=OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if OCMATFTE.objectivevaluecalc
        dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
    Jpar(:,OCMATFTE.optimalpartcoordinate(idx))=dxdt+Jt;
end
if OCMATFTE.continuationtype==1
    dxdt1=0;
    dxdt2=0;
    if any(arc==(OCMATFTE.continuationtimeindex-1))
        dxdt1=OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)*OCMATFTE.continuationvector(arc==(OCMATFTE.continuationtimeindex-1));
        if OCMATFTE.objectivevaluecalc
            dxdt1(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)*OCMATFTE.continuationvector(arc==(OCMATFTE.continuationtimeindex-1));
            %Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        end
        %Jpar(:,end)=dxdt+Jt;
    end
    if any(arc==OCMATFTE.continuationtimeindex)
        dxdt2=OCMATFTE.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)*OCMATFTE.continuationvector(arc==OCMATFTE.continuationtimeindex);
        if OCMATFTE.objectivevaluecalc
            dxdt2(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)*OCMATFTE.continuationvector(arc==OCMATFTE.continuationtimeindex);
            %Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        end
    end
        Jpar(:,end)=dxdt1-dxdt2;
end
if OCMATFTE.freeparameter
    Jmodelpar=dtds*OCMATFTE.canonicalsystemparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    Jpar(OCMATFTE.statecostatecoord{solutionindex},OCMATFTE.parametercoord)=Jmodelpar(:,OCMATFTE.parameterindex{solutionindex});
    if OCMATFTE.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jpar(OCMATFTE.objectivevaluecoord,OCMATFTE.parametercoord)=Jobjmodelpar(:,OCMATFTE.parameterindex{solutionindex});
    end
end
if OCMATFTE.continuationtype==2
    Jmodelpar=dtds*OCMATFTE.canonicalsystemparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);

    Jpar(OCMATFTE.statecostatecoord{solutionindex},end)=OCMATFTE.continuationvector(solutionindex)*Jmodelpar(:,OCMATFTE.parameterindex{solutionindex});
    if OCMATFTE.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATFTE.objectivefunctionparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jpar(OCMATFTE.objectivevaluecoord,end)=OCMATFTE.continuationvector(solutionindex)*Jobjmodelpar(:,OCMATFTE.parameterindex{solutionindex});
    end
end


%-------------------------------------------------------------------------
function res=bc2(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP

resconnec=[];
resinit=[];
restrans=[];
restarget=[];
respartconnect=[];

arctime=OCMATFTE.initialtimeinterval;
arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
end
if OCMATFTE.freeparameter
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=freepar(OCMATFTE.parametercoord);
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=OCMATFTE.initialparameter(ii,:)+freepar(end)*OCMATFTE.continuationvector(ii,:);
    end
end


if OCMATFTE.objectivevaluecalc
    OVal=0;
end
for ii=1:OCMATFTE.numberofmodels
    if ii==1
        initialstate=OCMATFTE.initialstate;
        if OCMATFTE.continuationtype==0
            initialstate=initialstate+freepar(end)*OCMATFTE.continuationvector;
        end
        restarget=depvara(OCMATFTE.statecoordinate{ii},1)-initialstate;
    end
    arcarg=OCMATCONT.HE.arcarg(OCMATFTE.arccoord{ii});
    if OCMATFTE.objectivevaluecalc
        OVal=OVal+OCMATFTE.salvagevalue{ii}(OCMATFTE.partendtime(ii),depvarb(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}(end)),modelpar{ii},arcarg(end));
    end

    for arc=1:OCMATFTE.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATFTE.reset{ii}(depvara(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),depvarb(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),modelpar{ii},arctime,arcarg,OCMATFTE.edge{ii},arc); ...
            OCMATFTE.guard{ii}(depvara(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),depvarb(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),modelpar{ii},arctime,arcarg,OCMATFTE.edge{ii},arc)];
        if OCMATFTE.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.objectivevaluecoord,OCMATFTE.arccoord{ii}(arc))-depvara(OCMATFTE.objectivevaluecoord,OCMATFTE.arccoord{ii}(arc+1))];
        end
    end
    if ii<OCMATFTE.numberofmodels
        respartconnect=[respartconnect; ...
            OCMATFTE.bcconnectingparts{ii}(depvara(OCMATFTE.statecostatecoord{ii},:),depvarb(OCMATFTE.statecostatecoord{ii},:),modelpar{ii},arctime,arcarg,[],OCMATFTE.cumsumnumarc(ii))];
        if OCMATFTE.optimalswitchingtime(ii)
            respartconnect=[respartconnect; ...
                OCMATFTE.bcoptimalconnectingparts{ii}(depvara(OCMATFTE.statecostatecoord{ii},:),depvarb(OCMATFTE.statecostatecoord{ii},:),modelpar,arctime,OCMATCONT.HE.arcarg,OCMATFTE.cumsumnumarc(ii))];
        end
        if OCMATFTE.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.objectivevaluecoord,OCMATFTE.cumsumnumarc(ii))-depvara(OCMATFTE.objectivevaluecoord,OCMATFTE.cumsumnumarc(ii)+1)];
        end
    end
end
if OCMATFTE.objectivevaluecalc
    resinit=[resinit; ...
        depvara(OCMATFTE.objectivevaluecoord,1)-OVal];
end
restrans=[restrans; ...
    OCMATFTE.bctransversality{OCMATFTE.numberofmodels}(OCMATFTE.partendtime(OCMATFTE.numberofmodels),depvarb(OCMATFTE.statecostatecoord{OCMATFTE.numberofmodels},OCMATFTE.arccoord{OCMATFTE.numberofmodels}(end)),modelpar{ii},arcarg(end))];

res=[resinit;restarget;restrans;resconnec;respartconnect];


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP

resconnec=[];
resinit=[];
restrans=[];
restarget=[];
respartconnect=[];

arctime=OCMATFTE.initialtimeinterval;
arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
end
if OCMATFTE.freeparameter
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=freepar(OCMATFTE.parametercoord);
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=OCMATFTE.initialparameter(ii,:)+freepar(end)*OCMATFTE.continuationvector(ii,:);
    end
end


if OCMATFTE.objectivevaluecalc
    OVal=0;
end
for ii=1:OCMATFTE.numberofmodels
    if ii==1
        initialstate=OCMATFTE.initialstate;
        if OCMATFTE.continuationtype==0
            initialstate=initialstate+freepar(end)*OCMATFTE.continuationvector;
        end
        restarget=depvara(OCMATFTE.statecoordinate{ii},1)-initialstate;
    end
    arcarg=OCMATCONT.HE.arcarg(OCMATFTE.arccoord{ii});
    if OCMATFTE.objectivevaluecalc
        OVal=OVal+OCMATFTE.salvagevalue{ii}(OCMATFTE.partendtime(ii),depvarb(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}(end)),modelpar{ii},arcarg(end));
    end

    for arc=1:OCMATFTE.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATFTE.reset{ii}(depvara(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),depvarb(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),modelpar{ii},arctime,arcarg,OCMATFTE.edge{ii},arc); ...
            OCMATFTE.guard{ii}(depvara(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),depvarb(OCMATFTE.statecostatecoord{ii},OCMATFTE.arccoord{ii}),modelpar{ii},arctime,arcarg,OCMATFTE.edge{ii},arc)];
        if OCMATFTE.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.objectivevaluecoord,OCMATFTE.arccoord{ii}(arc))-depvara(OCMATFTE.objectivevaluecoord,OCMATFTE.arccoord{ii}(arc+1))];
        end
    end
    if ii<OCMATFTE.numberofmodels
        respartconnect=[respartconnect; ...
            OCMATFTE.bcconnectingparts{ii}(depvara(OCMATFTE.statecostatecoord{ii},:),depvarb(OCMATFTE.statecostatecoord{ii},:),modelpar{ii},arctime,arcarg,OCMATFTE.edge{ii},OCMATFTE.cumsumnumarc(ii))];
        if OCMATFTE.optimalswitchingtime(ii)
            respartconnect=[respartconnect; ...
                OCMATFTE.bcoptimalconnectingparts{ii}(depvara(OCMATFTE.statecostatecoord{ii},:),depvarb(OCMATFTE.statecostatecoord{ii},:),modelpar,arctime,OCMATCONT.HE.arcarg,OCMATFTE.cumsumnumarc(ii))];
        end
        if OCMATFTE.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATFTE.objectivevaluecoord,OCMATFTE.cumsumnumarc(ii))-depvara(OCMATFTE.objectivevaluecoord,OCMATFTE.cumsumnumarc(ii)+1)];
        end
    end
end
if OCMATFTE.objectivevaluecalc
    resinit=[resinit; ...
        depvara(OCMATFTE.objectivevaluecoord,1)-OVal];
end
restrans=[restrans; ...
    OCMATFTE.bctransversality{OCMATFTE.numberofmodels}(OCMATFTE.partendtime(OCMATFTE.numberofmodels),depvarb(OCMATFTE.statecostatecoord{OCMATFTE.numberofmodels},OCMATFTE.arccoord{OCMATFTE.numberofmodels}(end)),modelpar{ii},arcarg(end))];

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

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATFTE.solutionindex(arc);
    arctime=OCMATFTE.initialtimeinterval;
    arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
    arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
    if OCMATFTE.continuationtype==1
        arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATFTE.testadmissibility{solutionindex}(t,sol.y(idx,:),modelpar{solutionindex},arcarg);
    else
        eqcoord=domainddata{solutionindex}(arcindex).eqcoord;
        [constr labelS]=OCMATFTE.testadmissibility{solutionindex}(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar{solutionindex},arcarg);
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
%figure(1)
h=OCMATFTE.plotcontinuation{1}(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE
idx=[];
if isempty(coeff)
    return
end
if OCMATFTE.findoptimalswitchingtime
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    arctime=OCMATFTE.initialtimeinterval;
    arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
    arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
    if OCMATFTE.continuationtype==1
        arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
    end
    solutionindex=OCMATFTE.solutionindex(OCMATFTE.findoptimalswitchingtime)-1;
    depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
    depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
    arcarg=OCMATCONT.HE.arcarg;
    out=OCMATFTE.bcoptimalconnectingparts{solutionindex}(depvara(OCMATFTE.statecostatecoord{solutionindex},:),depvarb(OCMATFTE.statecostatecoord{solutionindex},:),modelpar,arctime,arcarg,OCMATFTE.cumsumnumarc(solutionindex));
else
    out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
end
fprintf(1,' Distance to targer value parameter: %g\n',out);
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
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                arctime=OCMATFTE.initialtimeinterval;
                arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
                arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
                if OCMATFTE.continuationtype==1
                    arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
                end
                solutionindex=OCMATFTE.solutionindex(OCMATFTE.findoptimalswitchingtime)-1;
                depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                arcarg=OCMATCONT.HE.arcarg;
                out=OCMATFTE.bcoptimalconnectingparts{solutionindex}(depvara(OCMATFTE.statecostatecoord{solutionindex},:),depvarb(OCMATFTE.statecostatecoord{solutionindex},:),modelpar,arctime,arcarg,OCMATFTE.cumsumnumarc(solutionindex));
            elseif OCMATFTE.findoptimalparameter
                tpar=tangent(end);
                tangent=tangent/tpar;
                tangent=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
                out=tangent(OCMATFTE.objectivevaluecoord,end);
                OCMATFTE.derivative=out;

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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
if OCMATFTE.freeparameter
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=freepar(OCMATFTE.parametercoord);
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=OCMATFTE.initialparameter(ii,:)+freepar(end)*OCMATFTE.continuationvector(ii,:);
    end
end

arctime=OCMATFTE.initialtimeinterval;
arctime(OCMATFTE.optimalparttimeindex)=freepar(OCMATFTE.optimalpartcoordinate);
arctime(OCMATFTE.switchtimeindex)=freepar(OCMATFTE.switchtimecoordinate);
if OCMATFTE.continuationtype==1
    arctime(OCMATFTE.continuationtimeindex)=OCMATFTE.initialswitchtime+freepar(end)*OCMATFTE.continuationvector;
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
out.solverinfo.solutionindex=OCMATFTE.solutionindex;
out.solverinfo.initialstateindex=OCMATFTE.initialstateindex;
out.solverinfo.optimalswitchingtime=OCMATFTE.optimalswitchingtime;
out.solverinfo.continuationtype=OCMATFTE.continuationtype;
if OCMATFTE.continuationtype==0
elseif OCMATFTE.continuationtype==1
elseif OCMATFTE.continuationtype==2
    out.solverinfo.parameterindex=OCMATFTE.parameterindex;
    out.solverinfo.parameterindex=OCMATFTE.parameterindex;
    out.solverinfo.initialparameter=OCMATFTE.initialparameter;
    out.solverinfo.continuationvector=OCMATFTE.continuationvector;
end
 if OCMATFTE.objectivevaluecalc
     out.solverinfo.objectivevaluecoord=OCMATFTE.objectivevaluecoord;
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
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=freepar(OCMATFTE.parametercoord);
    end
end
if OCMATFTE.continuationtype==2
    for ii=1:OCMATFTE.numberofmodels
        modelpar{ii}(OCMATFTE.parameterindex{ii})=OCMATFTE.initialparameter(ii,:)+freepar(end)*OCMATFTE.continuationvector(ii,:);
    end
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


