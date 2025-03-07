function out=extremalt4ft()

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
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
end
arctime=OCMATFTE.initialarctimes;
arctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.freetimecoord).';
arctime(OCMATFTE.continuationtimeindex)=freepar(OCMATFTE.continuationtimecoordinate);
%arctime=[OCMATFTE.initialtime freepar(OCMATFTE.freetimecoord).'];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATFTE.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.vfreetimecoord);
    vdiffarctime=diff(varctime);

    vdts=vdiffarctime(arc);
    
    dxdt(OCMATFTE.variationaldynamicscoordinate,:)=dtds*OCMATFTE.variationaldynamics(t,depvar,modelpar,arcarg)+vdts*OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
end
if OCMATFTE.objectivevaluecalc
    dxdt(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
end
if OCMATFTE.exogenousfunction
    dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=dtds*OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
end
arctime=OCMATFTE.initialarctimes;
arctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.freetimecoord).'; 
arctime(OCMATFTE.continuationtimeindex)=freepar(OCMATFTE.continuationtimecoordinate);
%arctime=[OCMATFTE.initialtime freepar(OCMATFTE.freetimecoord).'];
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=[OCMATFTE.JX,OCMATFTE.Jext];
if OCMATFTE.variationalcalculation
    J0=OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
    J(OCMATFTE.dFDXcoord1,OCMATFTE.dFDXcoord2)=J0;
else
    J(OCMATFTE.dFDXcoord1,OCMATFTE.dFDXcoord2)=OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
end
if OCMATFTE.objectivevaluecalc
    J(OCMATFTE.dFODXcoord1,OCMATFTE.dFODXcoord2)=OCMATFTE.objectivefunctionjacobian(t,depvar,modelpar,arcarg);
end
if OCMATFTE.exogenousfunction
    J(OCMATFTE.dFEDXcoord1,OCMATFTE.dFEDXcoord2)=OCMATFTE.exogenousjacobian(t,depvar,modelpar,arcarg);
end
J=dtds*J;
if OCMATFTE.variationalcalculation
    varctime=zeros(1,length(arctime));
    varctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.vfreetimecoord);
    vdiffarctime=diff(varctime);

    vdts=vdiffarctime(arc);

    J(OCMATFTE.dFVDXcoord1,OCMATFTE.dFVDXcoord2)=dtds*OCMATFTE.variationaljacobian(t,depvar,modelpar,arcarg)+vdts*[J0 OCMATFTE.dFVDV];
end
Jpar=OCMATFTE.Jpar;
if OCMATCONT.HE.numarc>1
        dxdt=OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
        if ~OCMATFTE.autonomous
            Jt=OCMATFTE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
        else
            Jt=0;
        end
        if OCMATFTE.variationalcalculation
            dxdt(OCMATFTE.variationaldynamicscoordinate,:)=OCMATFTE.variationaldynamics(t,depvar,modelpar,arcarg);
            Jt(OCMATFTE.variationaldynamicscoordinate,:)=0;%OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
        end
        if OCMATFTE.objectivevaluecalc
            dxdt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
            Jt(OCMATFTE.objectivevaluecoord,:)=OCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
        end
        if OCMATFTE.exogenousfunction
            dxdt(OCMATFTE.exogenousdynamicscoordinate,:)=OCMATFTE.exogenousdynamics(t,depvar,modelpar,arcarg);
            Jt(OCMATFTE.exogenousdynamicscoordinate,:)=0;
        end
        idx=find(OCMATFTE.freetimeindex==arc);
        idxp1=find(OCMATFTE.freetimeindex==arc+1);
        if ~isempty(idx)
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.freetimecoord(idx))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
        if ~isempty(idxp1)
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.freetimecoord(idxp1))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        end
        idx=find(OCMATFTE.continuationtimeindex==arc);
        idxp1=find(OCMATFTE.continuationtimeindex==arc+1);
        if ~isempty(idx)
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.continuationtimecoordinate(idx))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
        if ~isempty(idxp1)
            Jpar(OCMATFTE.ODEcoord,OCMATFTE.continuationtimecoordinate(idxp1))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        end

else
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
    Jpar(:,OCMATFTE.continuationcoord)=dxdt+diffarctime(arc)*(s-arc)*Jt;
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
if ~isempty(OCMATFTE.freeparameter)
    modelpar(OCMATFTE.freeparameter)=freepar(OCMATFTE.freeparametercoordinate);
end
arctime=OCMATFTE.initialarctimes;
arctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.freetimecoord).'; 
arctime(OCMATFTE.continuationtimeindex)=freepar(OCMATFTE.continuationtimecoordinate);
switchtimes=arctime(2:end-1);
timehorizon=arctime(end);
resconnec=[];
userbc=[];

resinit=[OCMATFTE.bcinitial(depvara(:,1),OCMATFTE.fixinitstatecoord,OCMATFTE.initstate,modelpar,OCMATCONT.HE.arcarg(1));
    depvarb(OCMATFTE.fixendstatecoord,end)-OCMATFTE.endstate];
if OCMATFTE.variationalcalculation
    resinit=[resinit; ...
        OCMATFTE.variationalbcinitial(depvara,[],[],modelpar,OCMATCONT.HE.arcarg(1));];
end

if OCMATFTE.stateconstraint
    jumparg=freepar(OCMATFTE.entrytimecoordinate);
end
if OCMATFTE.variationalcalculation
    if OCMATFTE.stateconstraint
        vjumparg=freepar(OCMATFTE.ventrytimecoordinate);
    end
end
if OCMATFTE.userbc
    userbc=OCMATFTE.userbcfunc([OCMATFTE.initialtime switchtimes(:).' timehorizon],depvara,depvarb,modelpar,OCMATCONT.HE.edge);
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
    restrans(OCMATFTE.fixendstatecoord)=[];
end
if ~isempty(OCMATFTE.fixendcostatecoord)
    restrans(OCMATFTE.fixendcostatecoord)=depvarb(OCMATFTE.statenum+OCMATFTE.fixendcostatecoord,end)-OCMATFTE.endcostate;
else
    tmprestrans=OCMATFTE.bctransversality(0,depvara(:,1),modelpar,OCMATCONT.HE.arcarg(1));
    tmprestrans(OCMATFTE.fixinitstatecoord)=[];
    restrans=[restrans;tmprestrans];
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
    HT=OCMATFTE.startvalue+freepar(end)*OCMATFTE.continuationvector;
    restrans=[restrans(:); ...
        HT-OCMATFTE.bcoptimalhorizon(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end))];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    if OCMATFTE.stateconstraint && OCMATFTE.entryindex(ii)
        resconnec=[resconnec; ...
            OCMATFTE.bcstateconstraint(depvara,depvarb,modelpar,jumparg(OCMATFTE.entryindex(ii)),switchtimes,OCMATFTE.jumpid(ii),OCMATCONT.HE.edge,ii)];
    else
        resconnec=[resconnec;
            OCMATFTE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
    if any(OCMATFTE.freetimeindex==ii+1) %|| OCMATFTE.continuationtimeindex==ii+1
        resconnec=[resconnec; ...
            OCMATFTE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
    if OCMATFTE.variationalcalculation
        if OCMATFTE.stateconstraint && OCMATFTE.entryindex(ii)
            resconnec=[resconnec; ...
                OCMATFTE.variationalbcstateconstraint(depvara,depvarb,modelpar,[jumparg(OCMATFTE.entryindex(ii)) vjumparg(OCMATFTE.entryindex(ii))],switchtimes,OCMATFTE.jumpid(ii),OCMATCONT.HE.edge,ii)];
        else
            resconnec=[resconnec;
                OCMATFTE.variationalreset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
        end
        if any(OCMATFTE.freetimeindex==ii+1) || OCMATFTE.continuationtimeindex==ii+1
            resconnec=[resconnec; ...
                OCMATFTE.variationalguard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
        end
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
res=[resinit;resconnec;restrans(:);userbc];
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
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATFTE OCBVP
%domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    arctime=OCMATFTE.initialarctimes;
    arctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.freetimecoord).';
    arctime(OCMATFTE.continuationtimeindex)=freepar(OCMATFTE.continuationtimecoordinate);
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
arctime=OCMATFTE.initialarctimes;
arctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.freetimecoord).'; 
arctime(OCMATFTE.continuationtimeindex)=freepar(OCMATFTE.continuationtimecoordinate);
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
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
if OCMATFTE.findoptimalhorizon
    XT=y(:,end);
    T=freepar(OCMATFTE.continuationtimecoordinate(end));
    HT=OCMATFTE.bcoptimalhorizon(T,XT,modelpar,OCMATCONT.HE.arcarg(end));
    fprintf(1,' Hamiltonian at T: %g\n',HT);
end
if OCMATFTE.hitfunction
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
    depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);

    out=OCMATFTE.targetfunction(depvara,depvarb,modelpar,OCMATCONT.HE.arcarg);
    fprintf(1,' Hit value: %g\n',out);
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

failed=[];
for ii=id
    switch ii
        case 1
            if OCMATFTE.findoptimalhorizon
                XT=y(:,end);
                T=freepar(OCMATFTE.continuationtimecoordinate(end));
                out=OCMATFTE.bcoptimalhorizon(T,XT,modelpar,OCMATCONT.HE.arcarg(end));
            elseif OCMATFTE.hitfunction
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                out=OCMATFTE.targetfunction(depvara,depvarb,modelpar,OCMATCONT.HE.arcarg);
            else
                out=freepar(OCMATCONT.HE.numparameter)-OCMATFTE.targetvalue;
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
arctime=OCMATFTE.initialarctimes;
arctime(OCMATFTE.freetimeindex)=freepar(OCMATFTE.freetimecoord).'; 
arctime(OCMATFTE.continuationtimeindex)=freepar(OCMATFTE.continuationtimecoordinate);
%arctime=[OCMATFTE.initialtime freepar(OCMATFTE.freetimecoord).'];
out.arcinterval=arctime;
%out.arcinterval=[OCMATFTE.initialtime freepar(OCMATFTE.freetimecoord).' OCMATFTE.truncationtime];
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
out.solverinfo.conttype='extremalt4ft';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.freetimecoord=OCMATFTE.freetimecoord;
out.solverinfo.fixtimeindex=OCMATFTE.fixtimeindex;
out.solverinfo.objectivevaluecalc=OCMATFTE.objectivevaluecalc;
out.solverinfo.objectivevaluecoord=OCMATFTE.objectivevaluecoord;
out.solverinfo.fixendstatecoord=OCMATFTE.fixendstatecoord;
out.solverinfo.fixinitstatecoord=OCMATFTE.fixinitstatecoord;
if OCMATFTE.variationalcalculation
    out.solverinfo.vfreeparametercoordinate=OCMATFTE.vfreeparametercoordinate;
    out.solverinfo.ventrytimecoordinate=OCMATFTE.ventrytimecoordinate;
    out.solverinfo.vfreetimecoord=OCMATFTE.vfreetimecoord;
end

out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.solverinfo.stateconstraint=OCMATFTE.stateconstraint;
out.solverinfo.odenum=OCBVP.numode;

if OCMATFTE.stateconstraint
    out.solverinfo.entrytimecoordinate=OCMATFTE.entrytimecoordinate;
    out.solverinfo.entryindex=OCMATFTE.entryindex;
    out.solverinfo.jumpid=OCMATFTE.jumpid;
end
if OCMATFTE.exogenousfunction
    out.solverinfo.exogenousdynamicscoordinate=OCMATFTE.exogenousdynamicscoordinate;
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
z=[];
modelpar=OCMATFTE.parametervalue;
switch OCMATCONT.bvpmethod
    case 'gbvp4c'
        y=zeros(OCBVP.maxnumode,length(tmesh));
        y(OCMATCONT.HE.DDATA.meshvalcoord)=coeff(OCMATCONT.HE.ycoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
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
        save([OCMATFTE.basicglobalvarfilename '4extremalt4ft'],'MODELINFO')
    end
    save([OCMATFTE.basicresultfilename '4extremalt4ft'],'sout','bvpout')
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