function out=indifferencesolution4mm()

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
global OCMATINDIF


solutionindex=OCMATINDIF.solutionindex(arc);
if ~isempty(OCMATINDIF.freemodelparameterindex{solutionindex})
    modelpar{solutionindex}(OCMATINDIF.freemodelparameterindex{solutionindex})=freepar(OCMATINDIF.freemodelparametercoord{solutionindex});
end
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end

arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime(solutionindex)];

diffarctime=diff(arctime);
arcarg=OCMATINDIF.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
if OCMATINDIF.variationalcalculation(solutionindex)
    varctime=zeros(1,length(arctime));
    varctime(OCMATINDIF.freeswitchingtimeindex{solutionindex})=freepar(OCMATINDIF.vfreetimecoord{solutionindex});
    vdiffarctime=diff(varctime);
    vdts=vdiffarctime(relarc);
    dxdt(OCMATINDIF.variationaldynamicscoordinate{1},:)=dtds*OCMATINDIF.variationaldynamics{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)+vdts*OCMATINDIF.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
end

if OCMATINDIF.objectivevaluecalc(solutionindex)
    dxdt(OCMATINDIF.objectivevaluecoord(1),:)=dtds*OCMATINDIF.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if OCMATINDIF.includevariationalobjectivevalue(solutionindex)
        dxdt(OCMATINDIF.variationalobjectivevaluecoord(1),:)=dtds*OCMATINDIF.variationalobjectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)+vdts*OCMATINDIF.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT
solutionindex=OCMATINDIF.solutionindex(arc);
if ~isempty(OCMATINDIF.freemodelparameterindex{solutionindex})
    modelpar{solutionindex}(OCMATINDIF.freemodelparameterindex{solutionindex})=freepar(OCMATINDIF.freemodelparametercoord{solutionindex});
end
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end

arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime(solutionindex)];

diffarctime=diff(arctime);
arcarg=OCMATINDIF.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
J=OCMATINDIF.JX;
if OCMATINDIF.variationalcalculation(1)
    J0=OCMATINDIF.canonicalsystemjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg); % J0 is reused
    J(OCMATINDIF.dFDXcoord1,OCMATINDIF.dFDXcoord2)=dtds*J0;
else
    J(OCMATINDIF.dFDXcoord1,OCMATINDIF.dFDXcoord2)=dtds*OCMATINDIF.canonicalsystemjacobian{solutionindex}(t,depvar,modelpar,arcarg);
end
if OCMATINDIF.objectivevaluecalc(1)
    J(OCMATINDIF.dFODXcoord1(1),OCMATINDIF.dFODXcoord2)=dtds*OCMATINDIF.objectivefunctionjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if OCMATINDIF.includevariationalobjectivevalue(1)
        J(OCMATINDIF.dFODXcoord1(2),OCMATINDIF.dFODXcoord2)=dtds*OCMATINDIF.variationalobjectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
end
if OCMATINDIF.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(OCMATINDIF.freeswitchingtimeindex{solutionindex})=freepar(OCMATINDIF.vfreetimecoord{solutionindex});
    vdiffarctime=diff(varctime);
    vdts=vdiffarctime(relarc);
    %J(OCMATINDIF.dFVDXcoord1,OCMATINDIF.ODEcoord)=dtds*[OCMATINDIF.variationaljacobian(t,depvar,modelpar,arcarg) OCMATINDIF.dFDE]+vdts*[J0 OCMATINDIF.dFVDX OCMATINDIF.dFDO OCMATINDIF.dFDE];
    J(OCMATINDIF.dFVDXcoord1,OCMATINDIF.ODEcoord)=dtds*OCMATINDIF.variationaljacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)+vdts*[J0 OCMATINDIF.dFVDX OCMATINDIF.dFDO];
end
Jpar=OCMATINDIF.Jpar;
Jmodelpar=OCMATINDIF.canonicalsystemparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
if OCMATINDIF.numarc>1
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
        if OCMATINDIF.variationalcalculation
            if OCMATINDIF.includevariationalobjectivevalue
                dxdt(OCMATINDIF.variationalobjectivevaluecoord,:)=OCMATINDIF.variationalobjectivefunction(t,depvar,modelpar,arcarg);
                vdxdt(OCMATINDIF.objectivevaluecoord,:)=0;
                Jt(OCMATINDIF.variationalobjectivevaluecoord,:)=OCMATINDIF.variationalobjectivefunctionderivativetime(t,depvar,modelpar,arcarg);
                vdxdt(OCMATINDIF.variationalobjectivevaluecoord,:)=0;
            else
                vdxdt(OCMATINDIF.objectivevaluecoord,:)=0;
            end
        end
    end
    if OCMATINDIF.exogenousfunction
        if OCMATINDIF.variationalcalculation
            vdxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousvariationaldynamics(t,depvar,modelpar,arcarg);
        end
        dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
        Jt(OCMATINDIF.exogenousdynamicscoordinate,:)=0;
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if OCMATINDIF.variationalcalculation
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.vfreetimecoord(arc))=vdxdt;
        end
        if arc>1
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
            if OCMATINDIF.variationalcalculation
                Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.vfreetimecoord(arc-1))=-vdxdt;
            end
        end
    else
        Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        if OCMATINDIF.optimalhorizon
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        end
        if OCMATINDIF.variationalcalculation
            Jpar(OCMATINDIF.ODEcoord,OCMATINDIF.vfreetimecoord(arc-1))=-vdxdt;
        end
    end
else
end
if ~isempty(OCMATINDIF.freemodelparameterindex{solutionindex})
    Jpar(OCMATINDIF.statecostatecoord,OCMATINDIF.freemodelparametercoord{solutionindex})=dtds.*Jmodelpar(:,OCMATINDIF.freemodelparameterindex{solutionindex});
    if OCMATINDIF.variationalcalculation(1)
        Jvarmodelpar=dtds*OCMATINDIF.variationalparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jpar(OCMATINDIF.variationaldynamicscoordinate{1},OCMATINDIF.freemodelparametercoord{solutionindex})=Jvarmodelpar(:,OCMATINDIF.freemodelparameterindex{solutionindex})+vdts*Jmodelpar(:,OCMATINDIF.freemodelparameterindex{solutionindex});
    end
    if OCMATINDIF.objectivevaluecalc(1)
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        if OCMATINDIF.includevariationalobjectivevalue(1)
            Jvarobjmodelpar=dtds*OCMATINDIF.variationalobjectivefunctionparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        end
        Jpar(OCMATINDIF.objectivevaluecoord(1),OCMATINDIF.freemodelparametercoord{solutionindex})=Jobjmodelpar(:,OCMATINDIF.freemodelparameterindex{solutionindex});
        if OCMATINDIF.includevariationalobjectivevalue(1)
            Jpar(OCMATINDIF.variationalobjectivevaluecoord(1),OCMATINDIF.freemodelparametercoord{solutionindex})=Jvarobjmodelpar(:,OCMATINDIF.freemodelparameterindex{solutionindex});
        end
    end

end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF
resconnec=[];
resinit=[];
restrans=[];
resuser=[];
residpt=[];

O=zeros(1,OCMATINDIF.indifferenceorder);
freestatevector=0;
if ~isempty(OCMATINDIF.freestatevector)
    for jj=1:length(OCMATINDIF.freestatevectorcoord)
        freestatevector=freestatevector+freepar(OCMATINDIF.freestatevectorcoord(jj))*OCMATINDIF.freestatevector(:,jj);
    end
end
for ii=1:OCMATINDIF.indifferenceorder
    if ~isempty(OCMATINDIF.freemodelparameterindex{ii})
        modelpar{ii}(OCMATINDIF.freemodelparameterindex{ii})=freepar(OCMATINDIF.freemodelparametercoord{ii});
    end
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' OCMATINDIF.endtime(ii)];

    endtime=OCMATINDIF.endtime(ii);
    actdepvara=depvara(:,OCMATINDIF.arccoord{ii});
    actdepvarb=depvarb(:,OCMATINDIF.arccoord{ii});
    actarcarg=OCMATINDIF.arcarg(OCMATINDIF.arccoord{ii});
    actmodelpar=modelpar{ii};
    initialstate=OCMATINDIF.initialstate(:,ii);
    initialstate=initialstate+freestatevector;
    if ~isempty(OCMATINDIF.fixinitstate{ii})
        resinit=[resinit; ...
            actdepvara(OCMATINDIF.fixinitstate{ii},1)-OCMATINDIF.initialstate(OCMATINDIF.fixinitstate{ii},ii)];
        %initialstate(OCMATINDIF.fixinitstate{ii})=OCMATINDIF.initialstate(OCMATINDIF.fixinitstate{ii},ii);
    end
    if ii==1
        resinit=[resinit; ...
            actdepvara(OCMATINDIF.statecoord{1},1)-initialstate];
    end
    if OCMATINDIF.variationalcalculation(ii)
        resinit=[resinit; ...
            OCMATINDIF.variationalbcinitial{ii}(actdepvara(:,1),[],[],actmodelpar,actarcarg(1))];
    end
    if OCMATINDIF.variationalcalculation(ii)
        vswitchtimes=freepar(OCMATINDIF.vfreetimecoord{ii});
    end
    if OCMATINDIF.objectivevaluecalc(ii)
        OVal=OCMATINDIF.salvagevalue{ii}(endtime,actdepvarb(:,end),actmodelpar,actarcarg(end));
        resinit=[resinit; ...
            actdepvara(OCMATINDIF.objectivevaluecoord(1),1)-OVal];
        if OCMATINDIF.includevariationalobjectivevalue(ii)
            varOVal=OCMATINDIF.variationalsalvagevalue{ii}(endtime,actdepvarb(:,end),actmodelpar,actarcarg(end));
            resinit=[resinit;actdepvara(OCMATINDIF.variationalobjectivevaluecoord(1))-varOVal];
        end
    end
    if ii>1
        if ~isempty(OCMATINDIF.commoninitstate)
             resinit=[resinit; ...
                 actdepvara(OCMATINDIF.commoninitstate)-olddepvara(OCMATINDIF.commoninitstate)];
        end
    end
    restrans=[restrans; ...
        OCMATINDIF.bctransversality{ii}(endtime,actdepvarb([OCMATINDIF.statecoord{1} OCMATINDIF.costatecoord{1}],end),actmodelpar,actarcarg(end))];
    if OCMATINDIF.variationalcalculation(ii)
        restrans=[restrans; ...
            OCMATINDIF.variationalbctransversality{ii}(endtime,actdepvarb([OCMATINDIF.statecoord{1} OCMATINDIF.costatecoord{1} OCMATINDIF.variationaldynamicscoordinate{1}],end),actmodelpar,actarcarg(end));];
    end
    if OCMATINDIF.userbc(ii)
        resuser=[resuser; ...
            OCMATINDIF.userfunctionbc{ii}(arctime,actdepvara,actdepvarb,actmodelpar,actarcarg)];
    end
    for arc=1:OCMATINDIF.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset{ii}(actdepvara,actdepvarb,actmodelpar,freepar(OCMATINDIF.switchtimecoord{ii}).',actarcarg,OCMATINDIF.edge{ii},arc); ...
            OCMATINDIF.guard{ii}(actdepvara,actdepvarb,actmodelpar,freepar(OCMATINDIF.switchtimecoord{ii}).',actarcarg,OCMATINDIF.edge{ii},arc)];
        if OCMATINDIF.objectivevaluecalc(ii)
            resconnec=[resconnec; ...
                actdepvarb(OCMATINDIF.objectivevaluecoord(1),arc)-actdepvara(OCMATINDIF.objectivevaluecoord(1),arc+1)];
            if OCMATINDIF.includevariationalobjectivevalue(ii)
                resconnec=[resconnec; ...
                    actdepvarb(OCMATINDIF.variationalobjectivevaluecoord(1),arc)-actdepvara(OCMATINDIF.variationalobjectivevaluecoord(1),arc+1)];
            end
        end
        if OCMATINDIF.variationalcalculation(ii)
            resconnec=[resconnec; ...
                OCMATINDIF.variationalreset{ii}(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),actmodelpar,freepar(OCMATINDIF.switchtimecoord{ii}).',vswitchtimes,actarcarg,OCMATINDIF.edge{ii},arc); ...
                OCMATINDIF.variationalguard{ii}(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),actmodelpar,freepar(OCMATINDIF.switchtimecoord{ii}).',vswitchtimes,actarcarg,OCMATINDIF.edge{ii},arc)];
        end
    end
    if OCMATINDIF.objectivevaluecalc
        O(ii)=actdepvarb(OCMATINDIF.objectivevaluecoord(1),end);
    end
    olddepvara=actdepvara;
end
for ii=1:OCMATINDIF.indifferenceorder-1
        residpt=[residpt; ...
            O(ii)-O(ii+1)];
end

res=[resinit;restrans;resconnec;resuser;residpt];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATINDIF OCBVP

[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATINDIF.solutionindex(arc);
%     if ~isempty(OCMATINDIF.freemodelparameterindex{solutionindex})
%         modelpar{solutionindex}(OCMATINDIF.freemodelparameterindex{solutionindex})=freepar(OCMATINDIF.freemodelparametercoord{solutionindex});
%     end
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
    end
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime(solutionindex)];
    diffarctime=diff(arctime);
    arcarg=OCMATINDIF.arcarg(arc);
    if solutionindex>1
        transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        transformedtimeshift=0;
    end
    t=diffarctime(relarc)*(s(leftarcindex(arc):rightarcindex(arc))-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    [constr,labelS]=OCMATINDIF.testadmissibility{solutionindex}(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar{solutionindex},arcarg);
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
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
%sol=evalatmesh(tmesh,y,z);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
for arc=1:sum(OCMATINDIF.cumsumnumarc(end))
    solutionindex=OCMATINDIF.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
    end
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime(solutionindex)];
    diffarctime=diff(arctime);
    if solutionindex>1
        transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        transformedtimeshift=0;
    end
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(relarc)*(s(leftarcindex(arc):rightarcindex(arc))-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
end
h=OCMATINDIF.plotcontinuation{1}(t,y,modelpar{1},OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
    fprintf(1,' Difference to target value: %g\n',y(2,1)-OCMATINDIF.targetvalue);
    %fprintf(1,' Difference to target value: %g\n',modelpar{1}(9)-OCMATINDIF.targetvalue);



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
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                %out=modelpar{1}(9)-OCMATINDIF.targetvalue;
                out=y(2,1)-OCMATINDIF.targetvalue;
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
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' OCMATINDIF.endtime(ii)];
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
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
out.solverinfo.conttype='indifferencesolution4mm';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.indifferenceorder=OCMATINDIF.indifferenceorder;
out.solverinfo.statecostatecoord=OCMATINDIF.statecostatecoord;
out.solverinfo.objectivevaluecoord=OCMATINDIF.objectivevaluecoord;
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.userbc=OCMATINDIF.userbc;
out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out=rmfield(out,'yp');


%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATINDIF
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
y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
z=[];
freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
for ii=1:OCMATINDIF.indifferenceorder
    if ~isempty(OCMATINDIF.freemodelparameterindex{ii})
        modelpar{ii}(OCMATINDIF.freemodelparameterindex{ii})=freepar(OCMATINDIF.freemodelparametercoord{ii});
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
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4mm'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4mm'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATINDIF

pathname=OCMATINDIF.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;
