function out=indifferencesolution()

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
global OCMATINDIF OCMATCONT
if ~isempty(OCMATINDIF.freeparameter)
    modelpar(OCMATINDIF.freeparameter)=freepar(OCMATINDIF.freeparametercoordinate);
end
solutionindex=OCMATINDIF.solutionindex(arc);
limitcycleindex=OCMATINDIF.limitcycleindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if limitcycleindex
    arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate(limitcycleindex))];
else
    if OCMATINDIF.freeendtime(solutionindex)
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.truncationtimecoordinate(solutionindex))];
    else
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
    end
end
if solutionindex>1
    transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
else
    transformedtimeshift=0;
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
end
if OCMATINDIF.exogenousfunction
    dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=dtds*OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
solutionindex=OCMATINDIF.solutionindex(arc);
limitcycleindex=OCMATINDIF.limitcycleindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if limitcycleindex
    arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate(limitcycleindex))];
else
    if OCMATINDIF.freeendtime(solutionindex)
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.truncationtimecoordinate(solutionindex))];
    else
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));

J=OCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if length(OCBVP.numode)>1
    numode=OCBVP.numode(arc);
else
    numode=OCBVP.numode;
end
if OCMATINDIF.objectivevaluecalc
    J=[J; ...
        dtds*OCMATINDIF.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(numode,length(OCMATINDIF.objectivevaluecoord))];
end
Jpar=zeros(numode,OCMATCONT.HE.numparameter);
if OCMATINDIF.numarc(solutionindex)>1
    dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATINDIF.autonomous
        Jt=OCMATINDIF.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if OCMATINDIF.exogenousfunction
        dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
        if ~isempty(OCMATINDIF.exogenousderivativetime)
            Jet=OCMATINDIF.exogenousderivativetime(t,depvar,modelpar,arcarg);
            Jt(OCMATINDIF.exogenousdynamicscoordinate)=Jet;
        end
    end
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(1:numode,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc))=dxdt;
        if relarc>1
            Jpar(1:numode,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc-1))=-dxdt;
        end
    else
        Jpar(1:numode,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        if OCMATINDIF.freeendtime(solutionindex)
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.truncationtimecoordinate(solutionindex))=dxdt;
        end
    end
else
    if OCMATINDIF.freeendtime(solutionindex)
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        if OCMATINDIF.objectivevaluecalc
            dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
        end
        if OCMATINDIF.exogenousfunction
            dxdt(OCMATINDIF.exogenousdynamicscoordinate,:)=OCMATINDIF.exogenousdynamics(t,depvar,modelpar,arcarg);
        end
        Jpar(:,OCMATINDIF.truncationtimecoordinate(solutionindex))=dxdt;
    end
end
if ~isempty(OCMATINDIF.freeparameter)
    if OCMATINDIF.implicit
        Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg,OCMATINDIF.freeparameter);
    else
        Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
    end
    Jpar(:,OCMATINDIF.freeparametercoordinate)=Jmodelpar(:,OCMATINDIF.freeparameter);
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT
resinit=[];
resasym=[];
resconnec=[];
residpt=[];
restarget=[];
resricatti=[];
resequilibrium=[];
resper=[];
if ~isempty(OCMATINDIF.freeparameter)
    modelpar(OCMATINDIF.freeparameter)=freepar(OCMATINDIF.freeparametercoordinate);
end
if OCMATINDIF.limitcyclecounter
    % conditions for limitcycle
    initaldepvarlc=zeros(OCMATINDIF.statecostatecoordinate(end),OCMATINDIF.limitcyclecounter);
    for ii=1:OCMATINDIF.limitcyclecounter
        actdepvara=depvara(:,OCMATINDIF.limitcycleindex==ii);
        actdepvarb=depvarb(:,OCMATINDIF.limitcycleindex==ii);
        switchtimes=freepar(OCMATINDIF.switchtimecoordinate{OCMATINDIF.indifferenceorder+ii});
        arcarg=OCMATINDIF.arcarg{OCMATINDIF.indifferenceorder+ii};
        initaldepvarlc(:,ii)=actdepvara(:,1);
        resper=[resper; ...
            OCMATINDIF.bclimitcycle(actdepvara,actdepvarb); ...
            sum(OCMATINDIF.velocityvector{ii}.*(actdepvara(OCMATINDIF.velocitycoordinate{ii},1)-OCMATINDIF.initialpoint{ii}))];

        for arc=1:OCMATINDIF.numarc(OCMATINDIF.indifferenceorder+ii)-1
            resconnec=[resconnec; ...
                OCMATINDIF.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{OCMATINDIF.indifferenceorder+ii},arc); ...
                OCMATINDIF.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{OCMATINDIF.indifferenceorder+ii},arc)];
        end

    end
end
initialdepvar=zeros(size(depvara,1),OCMATINDIF.indifferenceorder);
if ~isempty(OCMATINDIF.player)
    enddepvar=zeros(size(depvara,1),OCMATINDIF.indifferenceorder);
end
initialarcarg=zeros(1,OCMATINDIF.indifferenceorder);
for ii=1:OCMATINDIF.indifferenceorder
    solutionindex=find(OCMATINDIF.solutionindex==ii);
    actdepvara=depvara(:,solutionindex);
    actdepvarb=depvarb(:,solutionindex);
    arcarg=OCMATINDIF.arcarg{ii};
    switchtimes=freepar(OCMATINDIF.switchtimecoordinate{ii});
    initialdepvar(:,ii)=actdepvara(:,1);
    if ~isempty(OCMATINDIF.player)
        enddepvar(:,ii)=actdepvarb(:,end);
    end
    initialarcarg(ii)=arcarg(1);
    if ii==1
        initialstate=OCMATINDIF.startvalue;
        if ~isempty(OCMATINDIF.freevectorcoordinate)
            for jj=1:length(OCMATINDIF.freevectorcoordinate)
                initialstate=initialstate+freepar(OCMATINDIF.freevectorcoordinate(jj))*OCMATINDIF.freevector(:,jj);
            end
            targetcoord=setdiff(OCMATINDIF.statecoordinate,OCMATINDIF.fixinitstate);

            %restarget=depvara(targetcoord,1)-initialstate(targetcoord);
            resinit=OCMATINDIF.bcinitial(depvara(:,1),targetcoord,initialstate(targetcoord),modelpar,arcarg(1));
        else
            if OCMATINDIF.implicit
                resinit=OCMATINDIF.algebraicequation(depvara(:,1),modelpar,arcarg(1));
            end

        end
    end
    if ii==1 && ~isempty(OCMATINDIF.fixinitstate)
        resinit=[resinit;...
            OCMATINDIF.bcinitial(initialdepvar(:,1),OCMATINDIF.fixinitstate,OCMATINDIF.startvalue(OCMATINDIF.fixinitstate),modelpar,arcarg(1))];
    end
    if ii>1
        try
            resinit=[resinit; ...
                OCMATINDIF.bcinitial(initialdepvar(:,ii),OCMATINDIF.statecoordinate,initialdepvar(OCMATINDIF.statecoordinate,ii-1),modelpar,arcarg(1))];
        catch
            resinit=[resinit; ...
                OCMATINDIF.bcinitial(initialdepvar(:,ii),OCMATINDIF.statecoordinate,initialdepvar(OCMATINDIF.statecoordinate,ii-1))];
        end
        if isempty(OCMATINDIF.player)
            residpt=[residpt; ...
                OCMATINDIF.bcindifference(initialdepvar,modelpar,initialarcarg,[ii-1 ii])];
        else
            residpt=[residpt; ...
                OCMATINDIF.bcindifference(enddepvar,modelpar,initialarcarg,[ii-1 ii],OCMATINDIF.objectivevaluecoord,OCMATINDIF.player)];
        end
    end
    if OCMATINDIF.objectivevaluecalc
        resinit=[resinit;actdepvara(OCMATINDIF.objectivevaluecoord)];
    end
    if OCMATINDIF.exogenousfunction
        resinit=[resinit; ...
            depvara(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(1))-OCMATINDIF.exogenousinitialstatesfunc(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,arcarg)];
        %         resinit=[resinit; ...
        %             depvara(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(1))-OCMATINDIF.exogenousinitialstates];
    end
    if ~isempty(OCMATINDIF.divergingcoordinate{ii})
        if ~isempty(OCMATINDIF.divergingcoordinate{ii})
            resasym=[resasym; ...
                OCMATINDIF.bcinf(depvarb(:,OCMATINDIF.cumsumnumarc(ii)),OCMATINDIF.fixendcoordinate{ii},OCMATINDIF.endvalue{ii},modelpar,arcarg(end))];
        end
    end
    if ~isempty(OCMATINDIF.freeparameter)
        if ~OCMATINDIF.simple(ii)
            hatx=freepar(OCMATINDIF.equilibriumcoordinate{ii});
            Y=freepar(OCMATINDIF.Ycoordinate{ii});
            asymptoticmatrix=OCMATINDIF.Q0{ii}*[-Y';OCMATINDIF.Id{ii}];
            Jac=OCMATINDIF.canonicalsystemjacobian(0,hatx,modelpar,OCMATINDIF.limsetarcarg{ii});
            if OCMATINDIF.implicit && OCMATINDIF.implicitcontrolnum(ii)
                dudx=OCMATINDIF.dimplicitcontroldx(0,hatx,modelpar,OCMATINDIF.limsetarcarg{ii});
                Jac=Jac(OCMATINDIF.statecostatecoordinate,OCMATINDIF.statecostatecoordinate)+Jac(OCMATINDIF.statecostatecoordinate,OCMATINDIF.statecostatecoordinate(end)+(1:OCMATINDIF.implicitcontrolnum(ii)))*dudx;
            end
            resricatti=[resricatti;ricatti(Y,Jac,OCMATINDIF.Q0{ii},OCMATINDIF.subspacedim{ii})];
        else
            hatx=freepar(OCMATINDIF.equilibriumcoordinate{ii});
            asymptoticmatrix=OCMATINDIF.asymptoticmatrix{ii};
        end
        if OCMATINDIF.equilibriumcounter>1 || ii==1
            resequilibrium=[resequilibrium;OCMATINDIF.equilibrium(hatx,modelpar,arcarg(end))];
        end

    else
        asymptoticmatrix=OCMATINDIF.asymptoticmatrix{ii};
        hatx=OCMATINDIF.saddlepoint{ii};
    end
    %try
    resasym=[resasym;OCMATINDIF.bcasymptotic(actdepvarb(:,end),asymptoticmatrix,hatx)];
    %catch

    %    resasym=[resasym;OCMATINDIF.bcasymptotic(depvarb(:,OCMATINDIF.arccoord{ii}),asymptoticmatrix,hatx,modelpar,arcarg(end))];
    %end
    if OCMATINDIF.fixdistance(ii)==1
        if isempty(OCMATINDIF.divergingcoordinate{ii})
            resasym=[resasym; ...
                sqrt(sum((OCMATINDIF.saddlepoint{ii}(OCMATINDIF.fixdistancecoordinate{ii})-depvarb(OCMATINDIF.fixdistancecoordinate{ii},OCMATINDIF.cumsumnumarc(ii))).^2))-OCMATINDIF.distance(ii)];
        else
            convergingcoordinate=setdiff(OCMATINDIF.statecostatecoordinate,OCMATINDIF.divergingcoordinate{ii});
            yend=depvarb(:,OCMATINDIF.arccoord{ii});
            hatx=OCMATINDIF.saddlepoint{ii};
            %yend(OCMATINDIF.divergingcoordinate{ii},:)=[];
            %hatx(OCMATINDIF.divergingcoordinate{ii},:)=[];
            resasym=[resasym; ...
                sqrt(sum((yend(convergingcoordinate)-hatx(convergingcoordinate)).^2))-OCMATINDIF.distance(ii)];
        end
    end
    for arc=1:OCMATINDIF.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
            OCMATINDIF.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
        if OCMATINDIF.objectivevaluecalc
            resconnec=[resconnec;actdepvara(OCMATINDIF.objectivevaluecoord,(arc+1))-actdepvarb(OCMATINDIF.objectivevaluecoord,(arc))];
        end
        if OCMATINDIF.exogenousfunction
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(arc))-depvara(OCMATINDIF.exogenousdynamicscoordinate,OCMATINDIF.arccoord{ii}(arc+1))];
        end
    end
end
res=[resinit;resasym;resconnec;residpt;resequilibrium;resricatti;resper];
%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATINDIF OCBVP

[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
if ~isempty(OCMATINDIF.freeparameter)
    modelpar(OCMATINDIF.freeparameter)=freepar(OCMATINDIF.freeparametercoordinate);
end

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATINDIF.solutionindex(arc);
    limitcycleindex=OCMATINDIF.limitcycleindex(arc);
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
    end
    if limitcycleindex
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate(limitcycleindex))];
    else
        if OCMATINDIF.freeendtime(solutionindex)
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.truncationtimecoordinate(solutionindex))];
        else
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
        end
    end
    diffarctime=diff(arctime);
    arcarg=OCMATCONT.HE.arcarg(arc);
    if solutionindex>1
        transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
    else
        transformedtimeshift=0;
    end
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr,labelS]=OCMATINDIF.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        %eqcoord=domainddata(arcindex).eqcoord;
        [constr,labelS]=OCMATINDIF.testadmissibility(t(leftarcindex(arc):rightarcindex(arc)),sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
    end
    violationmat=constr<-OCMATCONT.OPTIONS.admissibletol;
    if any(violationmat(:))
        counter=counter+1;
        [rows,cols]=find(violationmat);
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
for ii=1:OCMATINDIF.indifferenceorder
    checkcoordinate=setdiff(OCMATINDIF.statecostatecoordinate,OCMATINDIF.fixdistancecoordinate{ii});
    if ~isempty(checkcoordinate)
        if~isempty(OCMATINDIF.freeparameter)
            hatx=freepar(OCMATINDIF.equilibriumcoordinate{ii});
        else
            hatx=OCMATINDIF.saddlepoint{ii};
        end
        yend=y(OCMATINDIF.statecostatecoordinate,rightarcindex(OCMATINDIF.cumsumnumarc(ii)));
        if ~isempty(OCMATINDIF.divergingcoordinate)
            convergingcoordinate=setdiff(OCMATINDIF.statecostatecoordinate,OCMATINDIF.divergingcoordinate{ii});
            hatx=hatx(convergingcoordinate);
            yend=yend(convergingcoordinate);
        end
        violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend(checkcoordinate)-hatx(checkcoordinate))<0;
        if violationmat
            counter=counter+1;
            cols=size(y,2);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='maxdistance';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=norm(yend-hatx);
            infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx);
            b=min([b infoS(counter).minval]);
        end

    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
for arc=1:sum(OCMATINDIF.cumsumnumarc(end))
    solutionindex=OCMATINDIF.solutionindex(arc);
    limitcycleindex=OCMATINDIF.limitcycleindex(arc);
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
    end
    if limitcycleindex
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate(limitcycleindex))];
    else
        if OCMATINDIF.freeendtime(solutionindex)
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.truncationtimecoordinate(solutionindex))];
        else
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
        end
    end
    diffarctime=diff(arctime);
    if solutionindex>1
        transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
    else
        transformedtimeshift=0;
    end
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(relarc)*(s(leftarcindex(arc):rightarcindex(arc))-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
end
h=OCMATINDIF.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
if ~isempty(OCMATINDIF.targetparameter)
    fprintf(1,' Difference to target parameter value: %g\n',OCMATINDIF.targetparametervalue-modelpar(OCMATINDIF.targetparameter));
elseif OCMATINDIF.hitfunction
    depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
    depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
    fprintf(1,' Hit value: %g\n',OCMATINDIF.targetfunction(depvara,depvarb,modelpar,OCMATINDIF));
elseif ~isempty(OCMATINDIF.targetvalue)
    if ~isempty(OCMATINDIF.targetcoordinate)
        fprintf(1,' Difference to target value: %g\n',OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1));
    elseif ~isempty(OCMATINDIF.targetvectorcoordinate)
        fprintf(1,' Difference to target value: %g\n',OCMATINDIF.targetvalue-freepar(OCMATINDIF.freevectorcoordinate(OCMATINDIF.targetvectorcoordinate)));
    end
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
            [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
            if ~isempty(OCMATINDIF.targetparameter)
                out=OCMATINDIF.targetparametervalue-modelpar(OCMATINDIF.targetparameter);
            elseif OCMATINDIF.hitfunction
                depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                out=OCMATINDIF.targetfunction(depvara,depvarb,modelpar,OCMATINDIF);
            elseif~isempty(OCMATINDIF.targetcoordinate)
                out=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
            elseif ~isempty(OCMATINDIF.targetvectorcoordinate)
                out=OCMATINDIF.targetvalue-freepar(OCMATINDIF.freevectorcoordinate(OCMATINDIF.targetvectorcoordinate));
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
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[];
out.timehorizon=[];
for ii=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.freeendtime(ii)
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ii}).' freepar(OCMATINDIF.truncationtimecoordinate(ii))];
    else
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ii}).' OCMATINDIF.truncationtime(ii)];
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
    if ~isempty(OCMATINDIF.freeparameter)
        out.solverinfo.equilibriumcoordinate=OCMATINDIF.equilibriumcoordinate;
    if ~OCMATINDIF.simple(ii)
        out.solverinfo.Ycoordinate=OCMATINDIF.Ycoordinate;
        out.solverinfo.subspacedim=OCMATINDIF.subspacedim;
        out.solverinfo.orthspacedim=OCMATINDIF.orthspacedim;
        out.solverinfo.qbasis=OCMATINDIF.Q0;
    end
    end
    if OCMATINDIF.freeendtime(ii)
        out.solverinfo.truncationtimecoordinate=OCMATINDIF.truncationtimecoordinate;
    end
end
ctr=ii;
for ii=1:OCMATINDIF.limitcyclecounter
    ctr=ctr+1;
    arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ctr}).' freepar(OCMATINDIF.periodcoordinate(ii))];
    out.solverinfo.arcinterval{ctr}=arctime;
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATINDIF.initialtime;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATINDIF.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoordinate=OCMATINDIF.switchtimecoordinate;
out.solverinfo.divergingcoordinate=OCMATINDIF.divergingcoordinate;
%out.solverinfo.convergingcoordinate=OCMATINDIF.convergingcoordinate;
out.solverinfo.fixendcoordinate=OCMATINDIF.fixendcoordinate;
out.solverinfo.fixdistance=OCMATINDIF.fixdistance;
out.solverinfo.freeendtime=OCMATINDIF.freeendtime;
out.solverinfo.odenum=OCBVP.numode;
out.solverinfo.solutionindex=OCMATINDIF.solutionindex;
out.solverinfo.freevector=OCMATINDIF.freevector;
out.solverinfo.freevectorcoordinate=OCMATINDIF.freevectorcoordinate;
if ~isempty(OCMATINDIF.player)
    out.solverinfo.player=OCMATINDIF.player;
end

out.solverinfo.freeparameter=OCMATINDIF.freeparameter;
if OCMATINDIF.freeparameter
    out.solverinfo.freeparameter=OCMATINDIF.freeparameter;
    out.solverinfo.freeparametercoordinate=OCMATINDIF.freeparametercoordinate;
end

switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.solverinfo.yp=out.yp;
        out.solverinfo.ypmid=out.ypmid;
        out=rmfield(out,{'yp','ypmid'});
    case {'bvp4c','gbvp4'}
        out.solverinfo.yp=out.yp;
        out=rmfield(out,'yp');
end
if OCMATINDIF.stateconstraint
    out.solverinfo.jumpcostatecoord=OCMATINDIF.jumpcostatecoord;
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
    case 'gbvp4c'
        y=zeros(OCMATCONT.maxnumode,length(tmesh));
        y(OCMATCONT.HE.DDATA.meshvalcoord)=coeff(OCMATCONT.HE.ycoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        z=[];
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
if ~isempty(OCMATINDIF.freeparameter)
    modelpar(OCMATINDIF.freeparameter)=freepar(OCMATINDIF.freeparametercoordinate);
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
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution'],'sout','bvpout')
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
global OCMATCONT OCMATINDIF


if ~isempty(OCMATINDIF.freeparameter)
    for ii=1:OCMATINDIF.indifferenceorder
        if ~OCMATINDIF.simple(ii)
            [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
            Y=freepar(OCMATINDIF.Ycoordinate{ii});
            OCMATCONT.adapted = 1;
            % recalculate orthonormal basis and orthogonal complement
            [U,S,V]=svd(OCMATINDIF.Q0{ii}(:,1:OCMATINDIF.subspacedim{ii})+OCMATINDIF.Q0{ii}(:,OCMATINDIF.subspacedim{ii}+1:end)*Y);
            OCMATINDIF.Q0{ii}=U;
            OCMATINDIF.Y{ii}=zeros(OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});

            freepar(OCMATINDIF.Ycoordinate{ii})=OCMATINDIF.Y{ii};
            coeff=[y(:);freepar];
        else
            [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
            hatx=freepar(OCMATINDIF.equilibriumcoordinate{ii});
            Jac=OCMATINDIF.canonicalsystemjacobian(0,hatx,modelpar,OCMATINDIF.limsetarcarg{ii});
            if OCMATINDIF.implicit && OCMATINDIF.implicitcontrolnum(ii)
                dudx=OCMATINDIF.dimplicitcontroldx(0,hatx,modelpar,OCMATINDIF.limsetarcarg{ii});
                Jac=Jac(OCMATINDIF.statecostatecoordinate,OCMATINDIF.statecostatecoordinate)+Jac(OCMATINDIF.statecostatecoordinate,OCMATINDIF.statecostatecoordinate(end)+(1:OCMATINDIF.implicitcontrolnum(ii)))*dudx;
            end
            OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(Jac,OCMATINDIF.pathtype{ii},OCMATINDIF.maptype{ii},OCMATCONT.ZeroDeviationTolerance(ii),OCMATCONT.AsymptoticBCMethod{ii});
        end
    end
    flag = 1;
    coeff(OCMATCONT.HE.parametercoord)=freepar;
else
    flag=0;
end



%-----------------------------------------------------------------
function out=ricatti(Y,J,Q0,dim)
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=ricatticoefficient(Q0,J,dim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
