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
partindex=OCMATINDIF.arc2part(arc);
if OCMATINDIF.continuationtype==2
    modelpar{partindex}(OCMATINDIF.parameterindex{partindex})=OCMATINDIF.initialparameter{partindex}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{partindex};
end
if ~isempty(OCMATINDIF.freeparametervector)
    for jj=1:length(OCMATINDIF.freeparametercoordinate)
        modelpar{partindex}(OCMATINDIF.freeparametervectorcoordinate{jj}(:,partindex))=OCMATINDIF.freeinitialparametervalue{jj}(:,partindex)+freepar(OCMATINDIF.freeparametercoordinate(jj))*OCMATINDIF.freeparametervector{jj}(:,partindex);
    end
end
solutionindex=OCMATINDIF.arc2indifference(arc);
relarc=OCMATINDIF.arc2relarc(arc);
transformedtimeshift=OCMATINDIF.arc2reltimeshift(arc);
arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.numarcpersolution(solutionindex)+1:OCMATINDIF.numarcpersolution(solutionindex+1));
arctime(OCMATINDIF.switchtimeindex{solutionindex})=freepar(OCMATINDIF.switchtimecoordinate{solutionindex});
arctime(OCMATINDIF.optimalparttimeindex{solutionindex})=freepar(OCMATINDIF.opttimecoordinate{solutionindex});
if OCMATINDIF.continuationtype==1
    arctime(OCMATINDIF.continuationindex(1))=OCMATINDIF.initialarcinterval(OCMATINDIF.continuationindex(1))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(1);
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem{partindex}(t,depvar,modelpar{partindex},arcarg);
if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction{partindex}(t,depvar,modelpar{partindex},arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
solutionindex=OCMATINDIF.arc2indifference(arc);
partindex=OCMATINDIF.arc2part(arc);
relarc=OCMATINDIF.arc2relarc(arc);
transformedtimeshift=OCMATINDIF.arc2reltimeshift(arc);

arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.numarcpersolution(solutionindex)+1:OCMATINDIF.numarcpersolution(solutionindex+1));
arctime(OCMATINDIF.switchtimeindex{solutionindex})=freepar(OCMATINDIF.switchtimecoordinate{solutionindex});
arctime(OCMATINDIF.optimalparttimeindex{solutionindex})=freepar(OCMATINDIF.opttimecoordinate{solutionindex});
if OCMATINDIF.continuationtype==1
    arctime(OCMATINDIF.continuationindex(1))=OCMATINDIF.initialarcinterval(OCMATINDIF.continuationindex(1))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(1);
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);

t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATINDIF.canonicalsystemjacobian{partindex}(t,depvar,modelpar{partindex},arcarg);
J=dtds*J;
if OCMATINDIF.objectivevaluecalc
    J=[J; ...
        dtds*OCMATINDIF.objectivefunctionjacobian{partindex}(t,depvar,modelpar{partindex},arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA{partindex}(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA{partindex}(arcindex).numeq,OCMATCONT.HE.numparameter);
if any(arc==OCMATINDIF.switchtimeindex{solutionindex})
    dxdt=OCMATINDIF.canonicalsystem{partindex}(t,depvar,modelpar{partindex},arcarg);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{partindex}(t,depvar,modelpar{partindex},arcarg);
    end
    Jpar(:,OCMATINDIF.switchtimecoordinate(arc==OCMATINDIF.switchtimeindex))=-dxdt;
end
if any(arc==OCMATINDIF.switchtimeindex{solutionindex}-1)
    dxdt=OCMATINDIF.canonicalsystem{partindex}(t,depvar,modelpar{partindex},arcarg);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{partindex}(t,depvar,modelpar{partindex},arcarg);
    end
    Jpar(:,OCMATINDIF.switchtimecoordinate(arc==OCMATINDIF.switchtimeindex-1))=dxdt;
end
if OCMATINDIF.arcindex4optimaltime(arc)
    dxdt=OCMATINDIF.canonicalsystem{partindex}(t,depvar,modelpar{partindex},arcarg);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{partindex}(t,depvar,modelpar{partindex},arcarg);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{partindex}(t,depvar,modelpar{partindex},arcarg);
    end
    Jpar(:,OCMATINDIF.arc2optimaltimeindex(arc))=dxdt+Jt;
end
if arc>1 && OCMATINDIF.arcindex4optimaltime(arc-1)
    dxdt=OCMATINDIF.canonicalsystem{partindex}(t,depvar,modelpar{partindex},arcarg);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{partindex}(t,depvar,modelpar{partindex},arcarg);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{partindex}(t,depvar,modelpar{partindex},arcarg);
    end
    Jpar(:,OCMATINDIF.arc2optimaltimeindex(arc-1))=-dxdt+Jt;
end
if OCMATINDIF.continuationtype==2
    Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian{partindex}(t,depvar,modelpar{partindex},arcarg);
    Jpar(OCMATINDIF.statecostatecoord{partindex},OCMATINDIF.continuationcoordinate)=OCMATINDIF.continuationvector{partindex}*Jmodelpar(:,OCMATINDIF.parameterindex{partindex});
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian{partindex}(t,depvar,modelpar{partindex},arcarg);
        Jpar(OCMATINDIF.objectivevaluecoord,OCMATINDIF.continuationcoordinate)=OCMATINDIF.continuationvector{partindex}*Jobjmodelpar(:,OCMATINDIF.parameterindex{partindex});
    end
end
if OCMATINDIF.continuationtype==1 && any(arc==OCMATINDIF.conttimeindex2arc)
    idx=find(arc==OCMATINDIF.conttimeindex2arc);
    dxdt=OCMATINDIF.canonicalsystem{partindex}(t,depvar,modelpar{partindex},arcarg);
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{partindex}(t,depvar,modelpar{partindex},arcarg);
    end
    Jpar(:,OCMATINDIF.continuationcoordinate)=OCMATINDIF.continuationvector(idx)*dxdt;
end

if  ~isempty(OCMATINDIF.freeparametervector)
    Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian{partindex}(t,depvar,modelpar{partindex},arcarg);
    Jpar(OCMATINDIF.statecostatecoord{partindex},OCMATINDIF.freeparametercoordinate)=Jmodelpar(:,OCMATINDIF.freeparametervectorcoordinate{1}(arc));
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian{partindex}(t,depvar,modelpar{partindex},arcarg);
        Jpar(OCMATINDIF.objectivevaluecoord,OCMATINDIF.freeparametercoordinate)=Jobjmodelpar(:,OCMATINDIF.freeparametervectorcoordinate{1}(arc));
    end
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
    for ii=1:OCMATINDIF.numberofparts
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for jj=1:length(OCMATINDIF.freeparametercoordinate)
        for ii=1:OCMATINDIF.numberofparts
            modelpar{ii}(OCMATINDIF.freeparametervectorcoordinate{jj}(:,ii))=OCMATINDIF.freeinitialparametervalue{jj}(:,ii)+freepar(OCMATINDIF.freeparametercoordinate(jj))*OCMATINDIF.freeparametervector{jj}(:,ii);
        end
    end
end

O=zeros(1,OCMATINDIF.indifferenceorder);
for ii=1:OCMATINDIF.indifferenceorder
    if ii==1
        initialstate=OCMATINDIF.initialstate(:,1);
        if OCMATINDIF.continuationtype==0
            initialstate=initialstate+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector;
        end

        if ~isempty(OCMATINDIF.freestatevector)
            for jj=1:length(OCMATINDIF.freestatevectorcoordinate)
                initialstate=initialstate+freepar(OCMATINDIF.freestatevectorcoordinate(jj))*OCMATINDIF.freestatevector(:,jj);
            end
        end
        restarget=depvara(OCMATINDIF.statecoordinate{ii},1)-initialstate;
        if ~isempty(OCMATINDIF.fixinitialstate)
            restarget(OCMATINDIF.fixinitialstatecoordinate)=depvara(OCMATINDIF.fixinitialstatecoordinate,1)-OCMATINDIF.initialstate(OCMATINDIF.fixinitialstatecoordinate,1);
        end
    end

    depvarbcoord=OCMATINDIF.arcstructure{OCMATINDIF.partstructure{ii}(end)}(end);
    depvaracoord=OCMATINDIF.arcstructure{OCMATINDIF.partstructure{ii}(1)}(1);
    depvarpartcoord=OCMATINDIF.indiffence2partcoordinate{ii};
    if ii<OCMATINDIF.indifferenceorder
        depvaracoordp1=OCMATINDIF.arcstructure{OCMATINDIF.partstructure{ii+1}(1)}(1);
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate{ii},depvaracoord)-depvara(OCMATINDIF.statecoordinate{ii},depvaracoordp1)];
    end
    arcarg=OCMATINDIF.arcargument{OCMATINDIF.partstructure{ii}(end)};
    arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.numarcpersolution(ii)+1:OCMATINDIF.numarcpersolution(ii+1));
    if OCMATINDIF.continuationtype==1
        arctime(OCMATINDIF.continuationindex(1))=OCMATINDIF.initialarcinterval(OCMATINDIF.continuationindex(1))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(1);
    else
        arctime(OCMATINDIF.switchtimeindex{ii})=freepar(OCMATINDIF.switchtimecoordinate{ii});
        arctime(OCMATINDIF.optimalparttimeindex{ii})=freepar(OCMATINDIF.opttimecoordinate{ii});
    end
    switchtimes{ii}=arctime(2:end-1);
    restrans=[restrans; ...
        OCMATINDIF.bctransversality{OCMATINDIF.partstructure{ii}(end)}(arctime(end),depvarb(OCMATINDIF.statecostatecoord{ii},depvarbcoord),modelpar{OCMATINDIF.partstructure{ii}(end)},arcarg(end))];
    if OCMATINDIF.objectivevaluecalc
        OVal=OCMATINDIF.salvagevalue{OCMATINDIF.partstructure{ii}(end)}(arctime(end),depvarb(OCMATINDIF.statecostatecoord{ii},depvarbcoord),modelpar{OCMATINDIF.partstructure{ii}(end)},arcarg(end));
        resinit=[resinit; ...
            depvara(OCMATINDIF.objectivevaluecoord,depvaracoord)-OVal];
    end
    if OCMATINDIF.objectivevaluecalc
        O(ii)=depvarb(OCMATINDIF.objectivevaluecoord,depvarbcoord);
    end
    for jj=OCMATINDIF.partstructure{ii}(1:end-1)
        part=OCMATINDIF.parts2localpartsindex(jj);
        respartconnect=[respartconnect; ...
            OCMATINDIF.bcconnectingparts{jj}(depvara(OCMATINDIF.statecostatecoord{jj},depvarpartcoord),depvarb(OCMATINDIF.statecostatecoord{jj},depvarpartcoord),modelpar(OCMATINDIF.partstructure{ii}),arctime(OCMATINDIF.indiffence2parttimeindex{ii}),OCMATINDIF.indiffence2arcarguments{ii},part)];
        if OCMATINDIF.optimalswitchingtime(ii)
            respartconnect=[respartconnect; ...
                OCMATINDIF.bcoptimalconnectingparts{jj}(depvara(OCMATINDIF.statecostatecoord{jj},depvarpartcoord),depvarb(OCMATINDIF.statecostatecoord{jj},depvarpartcoord),modelpar(OCMATINDIF.partstructure{ii}),arctime(OCMATINDIF.indiffence2parttimeindex{ii}),OCMATINDIF.indiffence2arcarguments{ii},part)];
        end
        if OCMATINDIF.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.objectivevaluecoord,depvarpartcoord(part))-depvara(OCMATINDIF.objectivevaluecoord,depvarpartcoord(part+1))];
        end
    end
end
for ii=1:OCMATINDIF.numberofparts
    solutionindex=OCMATINDIF.arc2indifference(ii);
    arcarg=OCMATINDIF.arcargument{ii};

    for arc=1:OCMATINDIF.parts2numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset{ii}(depvara(:,OCMATINDIF.parts2arccoordinate{ii}),depvarb(:,OCMATINDIF.parts2arccoordinate{ii}),modelpar{ii},switchtimes{solutionindex}(OCMATINDIF.parts2switchtimeindex{ii}),arcarg,OCMATINDIF.edge{ii},arc); ...
            OCMATINDIF.guard{ii}(depvara(:,OCMATINDIF.parts2arccoordinate{ii}),depvarb(:,OCMATINDIF.parts2arccoordinate{ii}),modelpar{ii},switchtimes{solutionindex},arcarg,OCMATINDIF.edge{ii},arc)];
        if OCMATINDIF.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.objectivevaluecoord,OCMATINDIF.parts2arccoordinate{ii}(arc))-depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.parts2arccoordinate{ii}(arc+1))];
        end
    end
end
for ii=1:OCMATINDIF.indifferenceorder-1
    residpt=[residpt; ...
        O(ii)-O(ii+1)];
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
    solutionindex=OCMATINDIF.arc2indifference(arc);
    relarc=OCMATINDIF.arc2relarc(arc);
    transformedtimeshift=OCMATINDIF.arc2reltimeshift(arc);
    arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.numarcpersolution(solutionindex)+1:OCMATINDIF.numarcpersolution(solutionindex+1));
    if OCMATINDIF.continuationtype==1
        arctime(OCMATINDIF.continuationindex(1))=OCMATINDIF.initialarcinterval(OCMATINDIF.continuationindex(1))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(1);
    else
        arctime(OCMATINDIF.switchtimeindex{solutionindex})=freepar(OCMATINDIF.switchtimecoordinate{solutionindex});
        arctime(OCMATINDIF.optimalparttimeindex{solutionindex})=freepar(OCMATINDIF.opttimecoordinate{solutionindex});
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATINDIF.testadmissibility{solutionindex}(t,sol.y(idx,:),modelpar{solutionindex},arcarg);
    else
        eqcoord=domainddata{solutionindex}(arcindex).eqcoord;
        [constr labelS]=OCMATINDIF.testadmissibility{solutionindex}(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar{solutionindex},arcarg);
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
h=OCMATINDIF.plotcontinuation{1}(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
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
global OCMATCONT OCMATINDIF

failed=[];
out=[];
for ii=id
    switch ii
        case 1
            out=1-coeff(end);
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
for ii=1:OCMATINDIF.indifferenceorder
    arctime=OCMATINDIF.initialarcinterval(OCMATINDIF.numarcpersolution(ii)+1:OCMATINDIF.numarcpersolution(ii+1));
    arctime(OCMATINDIF.switchtimeindex{ii})=freepar(OCMATINDIF.switchtimecoordinate{ii});
    arctime(OCMATINDIF.optimalparttimeindex{ii})=freepar(OCMATINDIF.opttimecoordinate{ii});
    if OCMATINDIF.continuationtype==1
        arctime(OCMATINDIF.continuationindex(1))=OCMATINDIF.initialarcinterval(OCMATINDIF.continuationindex(1))+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector(1);
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
end
if OCMATINDIF.continuationtype==2
    for ii=1:OCMATINDIF.numberofparts
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for jj=1:length(OCMATINDIF.freeparametercoordinate)
        for ii=1:OCMATINDIF.numberofparts
            modelpar{ii}(OCMATINDIF.freeparametervectorcoordinate{jj}(:,ii))=OCMATINDIF.freeinitialparametervalue{jj}(:,ii)+freepar(OCMATINDIF.freeparametercoordinate(jj))*OCMATINDIF.freeparametervector{jj}(:,ii);
        end
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
out.solverinfo.initcoord=OCMATINDIF.initcoord;
out.solverinfo.solutionindex=OCMATINDIF.arc2indifference;
out.solverinfo.indifferenceorder=OCMATINDIF.indifferenceorder;
out.solverinfo.switchtimecoordinate=OCMATINDIF.switchtimecoordinate;
out.solverinfo.statecostatecoord=OCMATINDIF.statecostatecoord;
out.solverinfo.partindex=OCMATINDIF.partindex;
out.solverinfo.partstructure=OCMATINDIF.partstructure;
if OCMATINDIF.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATINDIF.objectivevaluecoord;
end
if OCMATINDIF.continuationtype==2
    out.solverinfo.parameterindex=OCMATINDIF.parameterindex;
end
out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;

out.solverinfo.numarc=OCMATINDIF.numarc;
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
    for ii=1:OCMATINDIF.numberofparts
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=OCMATINDIF.initialparameter{ii}+freepar(OCMATINDIF.continuationcoordinate)*OCMATINDIF.continuationvector{ii};
    end
end
if ~isempty(OCMATINDIF.freeparametervector)
    for jj=1:length(OCMATINDIF.freeparametercoordinate)
        for ii=1:OCMATINDIF.numberofparts
            modelpar{ii}(OCMATINDIF.freeparametervectorcoordinate{jj}(:,ii))=OCMATINDIF.freeinitialparametervalue{jj}(:,ii)+freepar(OCMATINDIF.freeparametercoordinate(jj))*OCMATINDIF.freeparametervector{jj}(:,ii);
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