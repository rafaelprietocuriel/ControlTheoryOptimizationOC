function out=heteroclinicep2lc()

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
%calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
%[J,M]=calc_RHSJac4LC(t,y,z,freepar,modelpar,odefun,bcfun,icfunc,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
%OCMATCONT.monodromy=M;

res=calc_RHS(t,y,z,freepar,modelpar,odefun,bcfun,icfunc);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
[J,M]=calc_RHSJac4LC(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
OCMATCONT.monodromy=M;

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
global OCMATHET
solutionindex=OCMATHET.solutionindex(arc);
limitcycleindex=OCMATHET.limitcycleindex(arc);
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if limitcycleindex
    arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(limitcycleindex))];
else
    if OCMATHET.freeendtime(solutionindex)
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
    else
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATHET.arcarg{solutionindex}(relarc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
solutionindex=OCMATHET.solutionindex(arc);
limitcycleindex=OCMATHET.limitcycleindex(arc);
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if limitcycleindex
    arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(limitcycleindex))];
else
    if OCMATHET.freeendtime(solutionindex)
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
    else
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATHET.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if length(OCBVP.numode)>1
    numode=OCBVP.numode(arc);
else
    numode=OCBVP.numode;
end
Jpar=zeros(numode,OCMATCONT.HE.numparameter);
if OCMATHET.numarc(solutionindex)>1
    dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATHET.autonomous
        Jt=OCMATHET.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if relarc<OCMATHET.numarc(solutionindex)
        Jpar(1:numode,OCMATHET.switchtimecoordinate{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-relarc)*Jt;
        if relarc>1
            Jpar(1:numode,OCMATHET.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        end
    else
        Jpar(1:numode,OCMATHET.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        if limitcycleindex
            dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATHET.periodcoordinate(limitcycleindex))=dxdt;
        else
            if OCMATHET.freeendtime(solutionindex)
                dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
                Jpar(:,OCMATHET.truncationtimecoordinate(solutionindex))=dxdt;
            end
        end
    end
else
    if limitcycleindex
        dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATHET.periodcoordinate(limitcycleindex))=dxdt;
    else
        if OCMATHET.freeendtime(solutionindex)
            dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATHET.truncationtimecoordinate(solutionindex))=dxdt;
        end
    end
end
if OCMATHET.implicit
    Jmodelpar=dtds*OCMATHET.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg,OCMATHET.parameterindex);
else
    Jmodelpar=dtds*OCMATHET.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
end
Jpar(:,OCMATHET.parametervaluecoord)=Jmodelpar(:,OCMATHET.parameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT
resasym=[];
resinit=[];
resconnec=[];
resequilibrium=[];
resricatti=[];
resper=[];
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

% conditions for limitcycle
initaldepvarlc=zeros(OCMATHET.statecostatecoordinate(end),OCMATHET.limitcyclecounter);
for ii=1:OCMATHET.limitcyclecounter
    actdepvara=depvara(:,OCMATHET.limitcycleindex==ii);
    actdepvarb=depvarb(:,OCMATHET.limitcycleindex==ii);
    switchtimes=freepar(OCMATHET.switchtimecoordinate{OCMATHET.hetorder+ii});
    arcarg=OCMATHET.arcarg{OCMATHET.hetorder+ii};
    initaldepvarlc(:,ii)=actdepvara(:,1);
    resper=[resper; ...
        OCMATHET.bclimitcycle(actdepvara,actdepvarb)];
    if ~OCMATHET.fixlccoordinate(ii)
        resper=[resper; ...
            sum(OCMATHET.velocityvector{ii}.*(actdepvara(OCMATHET.velocitycoordinate{ii},1)-OCMATHET.initialpoint{ii}))];
    else
        resper=[resper; ...
            actdepvara(OCMATHET.fixlccoordinate(ii),1)-OCMATHET.fixlcvalue(ii)];
    end

    
    for arc=1:OCMATHET.numarc(OCMATHET.hetorder+ii)-1
        resconnec=[resconnec; ...
            OCMATHET.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATHET.edge{OCMATHET.hetorder+ii},arc); ...
            OCMATHET.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATHET.edge{OCMATHET.hetorder+ii},arc)];
    end
    
end
initialdepvar=zeros(size(depvara,1),OCMATHET.hetorder);
ctr=0;
for ii=1:OCMATHET.hetorder
    solutionindex=find(OCMATHET.solutionindex==ii);
    actdepvara=depvara(:,solutionindex);
    actdepvarb=depvarb(:,solutionindex);
    arcarg=OCMATHET.arcarg{ii};
    switchtimes=freepar(OCMATHET.switchtimecoordinate{ii});

    initialdepvar(:,ii)=actdepvara(:,1);
    switch OCMATHET.limitsettype{OCMATHET.octrajectory2limset(ii,2)}
        case 'e'
            if ~OCMATHET.simple(ii)
                hatx=freepar(OCMATHET.equilibriumcoordinate{OCMATHET.octrajectory2limset(ii,3)});
                Y=freepar(OCMATHET.Ycoordinate{ii});
                asymptoticmatrix=OCMATHET.Q0{ii}*[-Y';OCMATHET.Id{ii}];
                Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,OCMATHET.limsetarcarg{ii});
                if OCMATHET.implicit && OCMATHET.implicitcontrolnum(ii)
                    dudx=OCMATHET.dimplicitcontroldx(0,hatx,modelpar,OCMATHET.limsetarcarg{ii});
                    Jac=Jac(OCMATHET.statecostatecoordinate,OCMATHET.statecostatecoordinate)+Jac(OCMATHET.statecostatecoordinate,OCMATHET.statecostatecoordinate(end)+(1:OCMATHET.implicitcontrolnum(ii)))*dudx;
                end
                resricatti=[resricatti;ricatti(Y,Jac,OCMATHET.Q0{ii},OCMATHET.subspacedim{ii})];
            else
                hatx=freepar(OCMATHET.equilibriumcoordinate{OCMATHET.octrajectory2limset(ii,3)});
                asymptoticmatrix=OCMATHET.asymptoticmatrix{ii};
            end
            if OCMATHET.equilibriumcounter>1 || ii==1
                resequilibrium=[resequilibrium;OCMATHET.equilibrium(hatx,modelpar,arcarg(end))];
            end
        case 'l'
            hatx=initaldepvarlc(:,OCMATHET.octrajectory2limset(ii,3));
            if ~OCMATHET.simple(ii)
                Y=freepar(OCMATHET.Ycoordinate{ii});
                asymptoticmatrix=OCMATHET.Q0{ii}*[-Y';OCMATHET.Id{ii}];
                Jac=OCMATCONT.monodromy;%{OCMATHET.octrajectory2limset(ii,3)};
                if OCMATHET.implicit && OCMATHET.implicitcontrolnum(ii)
                    dudx=OCMATHET.dimplicitcontroldx(0,hatx,modelpar,OCMATHET.limsetarcarg{ii});
                    Jac=Jac(OCMATHET.statecostatecoordinate,OCMATHET.statecostatecoordinate)+Jac(OCMATHET.statecostatecoordinate,OCMATHET.statecostatecoordinate(end)+(1:OCMATHET.implicitcontrolnum(ii)))*dudx;
                end
                resricatti=[resricatti;ricatti(Y,Jac,OCMATHET.Q0{ii},OCMATHET.subspacedim{ii})];
            else
                asymptoticmatrix=OCMATHET.asymptoticmatrix{ii};
            end
    end
    resasym=[resasym;OCMATHET.bcasymptotic(actdepvarb(:,end),asymptoticmatrix,hatx)];
    if OCMATHET.fixdistance(ii)
        resasym=[resasym; ...
            sqrt(sum((hatx(OCMATHET.statecostatecoordinate)-actdepvarb(OCMATHET.statecostatecoordinate,end)).^2))-OCMATHET.distance(ii)];
    end
    for arc=1:OCMATHET.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATHET.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc); ...
            OCMATHET.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc)];
    end
    if mod(ii,2)==0
        ctr=ctr+1;
        % resinit=[resinit; ...
        %     actdepvara(OCMATHET.fixcoordinate,1)-OCMATHET.fixvalue(:,ctr)];

        if OCMATHET.findinitconnection
            initialdepvar(:,2)=actdepvara(:,1);
            % resinit=[resinit; ...
            %     initialdepvar(OCMATHET.statecostatecoordinate,2)-initialdepvar(OCMATHET.statecostatecoordinate,1)-freepar(end)*OCMATHET.initialstatedifference(:,ctr)];
            resinit=[resinit; ...
                initialdepvar(OCMATHET.freestatecoordinate,ii)-initialdepvar(OCMATHET.freestatecoordinate,ii-1); ...
                initialdepvar(OCMATHET.costatecoordinate,2)-initialdepvar(OCMATHET.costatecoordinate,1)-freepar(end)*OCMATHET.initialcostatedifference(:,ctr)];
        else
            initialdepvar(OCMATHET.statecostatecoordinate,2)=actdepvara(OCMATHET.statecostatecoordinate,1);
            resinit=[resinit; ...
                initialdepvar(OCMATHET.freestatecoordinate,ii)-initialdepvar(OCMATHET.freestatecoordinate,ii-1); ...
                initialdepvar(OCMATHET.costatecoordinate,2)-initialdepvar(OCMATHET.costatecoordinate,1)];
        end
    else
        initialdepvar(:,1)=actdepvara(:,1);
    end
    resinit=[resinit; ...
        OCMATHET.bcinitial(initialdepvar(:,ii),OCMATHET.fixcoordinate,OCMATHET.fixvalue,modelpar,arcarg(1))];
end
res=[resinit;resasym;resconnec;resequilibrium;resricatti;resper];
%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATHET OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
depvara=y(:,leftarcindex);
depvarb=y(:,rightarcindex);
depvaralc=depvara(:,OCMATHET.limitcycleindex==1);

counter=0;
ctr=1;

for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATHET.solutionindex(arc);
    limitcycleindex=OCMATHET.limitcycleindex(arc);
    if solutionindex>1
        relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    if limitcycleindex
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(limitcycleindex))];
    else
        if arc>OCMATHET.cumsumnumarc(solutionindex)
            ctr=ctr+1;
        end
        if OCMATHET.freeendtime(solutionindex)
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
        else
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
        end
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if 0%~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr,labelS]=OCMATHET.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        [constr,labelS]=OCMATHET.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    reverseflag=0;
    if ~limitcycleindex
        if strcmp(OCMATHET.pathtype{OCMATHET.solutionindex(arc)},'u')
            reverseflag=1;
        end
    end
    if ~reverseflag
        violationmat=diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    else
        violationmat=-diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    end
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
for ii=1:OCMATHET.hetorder
    %actdepvara=depvara(:,OCMATHET.solutionindex==ii);
    actdepvarb=depvarb(:,OCMATHET.solutionindex==ii);
    if ~OCMATHET.freeendtime(ii)
        switch OCMATHET.limitsettype{OCMATHET.octrajectory2limset(ii,2)}
            case 'e'
                hatx=freepar(OCMATHET.equilibriumcoordinate{OCMATHET.octrajectory2limset(ii,3)});
                hatx=hatx(OCMATHET.statecostatecoordinate);
            case 'l'
                hatx=depvaralc(:,1);
        end
        violationmat=OCMATCONT.OPTIONS.maxdistance-norm(actdepvarb(OCMATHET.statecostatecoordinate,end)-hatx)<0;
        if violationmat
            counter=counter+1;
            cols=size(y,2);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='maxdistance';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=norm(y(OCMATHET.statecostatecoordinate,rightarcindex(OCMATHET.cumsumnumarc(ii)))-hatx);
            infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(y(OCMATHET.statecostatecoordinate,rightarcindex(OCMATHET.cumsumnumarc(ii)))-hatx);
            b=min([b infoS(counter).minval]);
        end
    end
%     if OCMATHET.hopf(ii)
%         hatx=freepar(OCMATHET.equilibriumcoordinate{OCMATHET.octrajectory2limset(ii,3)});
%         Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,arcarg);
%         violationmat=any(abs(imag(eig(Jac)))>0);
%         if violationmat
%             imJac=max(abs(imag(eig(Jac))));
%             counter=counter+1;
%             cols=size(y,2);
%             infoS(counter).arcarg=arcarg;
%             infoS(counter).arcnum=arc;
%             infoS(counter).rows='hopf';
%             infoS(counter).cols=cols;
%             infoS(counter).violationmat=violationmat;
%             infoS(counter).constraintvalue=imJac;
%             infoS(counter).minval=-imJac;
%             b=min([b infoS(counter).minval]);
%         end
%     end
end

%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
ctr=0;
clf
t=[];
for ii=1:OCMATHET.hetorder
    solutionindex=find(OCMATHET.solutionindex==ii);
    limitcycleindex=OCMATHET.limitcycleindex(solutionindex(1));
    if ~limitcycleindex
        ctr=ctr+1;
        transformedtimeshift=OCMATHET.cumsumnumarc(ctr);
        if OCMATHET.freeendtime(ii)
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{ii}).' freepar(OCMATHET.truncationtimecoordinate(ii))];
        else
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{ii}).' OCMATHET.truncationtime(ii)];
        end
        diffarctime=diff(arctime);
        for arc=1:OCMATHET.numarc(ii)
            t=[t diffarctime(arc)*(s(leftarcindex(solutionindex(arc)):rightarcindex(solutionindex(arc)))-transformedtimeshift)+(arctime(arc)-diffarctime(arc)*(arc-1))];
        end
    end
end
h=OCMATHET.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
idx=[];
if isempty(coeff)
    return
end
if isempty(OCMATHET.targetparametervalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(end));
else
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
    fprintf(1,' Difference to targetvalue: %g\n',OCMATHET.targetparametervalue-modelpar(OCMATHET.targetparameterindex));
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET

failed=[];
for ii=id
    switch ii
        case 1
                if isempty(OCMATHET.targetparametervalue)
                    out=coeff(end);
                else
                    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                    out=OCMATHET.targetparametervalue-modelpar(OCMATHET.targetparameterindex);
                end
            %out=tangent(end);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATHET OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
out.x=t;
out.y=y;
out.arcinterval=[];
out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
    OCMATCONT.HE.TIMEDDATA.rightarcindex];

for ii=1:OCMATHET.hetorder
    if OCMATHET.freeendtime(ii)
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{ii}).' freepar(OCMATHET.truncationtimecoordinate(ii))];
    else
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{ii}).' OCMATHET.truncationtime(ii)];
    end
    out.solverinfo.arcinterval{ii}=arctime;
    out.solverinfo.timehorizon(ii)=arctime(end);
    if ~OCMATHET.simple(ii)
        out.solverinfo.Ycoordinate=OCMATHET.Ycoordinate;
        out.solverinfo.subspacedim=OCMATHET.subspacedim;
        out.solverinfo.orthspacedim=OCMATHET.orthspacedim;
        out.solverinfo.qbasis=OCMATHET.Q0;
    end
    if OCMATHET.freeendtime(ii)
        out.solverinfo.truncationtimecoordinate=OCMATHET.truncationtimecoordinate;
    end
end
ctr=ii;
for ii=1:OCMATHET.limitcyclecounter
    ctr=ctr+1;
    arctime=[0 freepar(OCMATHET.switchtimecoordinate{ctr}).' freepar(OCMATHET.periodcoordinate(ii))];
    out.solverinfo.arcinterval{ctr}=arctime;
end
out.arcinterval=arctime;
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATHET.initialtime;
%out.modelparameter=OCMATHET.parametervalue;
%out.modelparameter(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='heteroclinicep2lc';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATHET.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.arcarg=OCMATHET.arcarg;
out.solverinfo.switchtimecoordinate=OCMATHET.switchtimecoordinate;
out.solverinfo.implicit=OCMATHET.implicit;
out.solverinfo.freeparameter=OCMATHET.parameterindex;
out.solverinfo.freeparametercoordinate=OCMATHET.parametervaluecoord;
if OCMATHET.implicit
    out.solverinfo.odenum=OCMATCONT.numode;
    out.solverinfo.maxnumode=OCMATCONT.maxnumode;
end
out.solverinfo.equilibriumcoordinate=OCMATHET.equilibriumcoordinate;
if OCMATHET.limitcyclecounter
    out.solverinfo.periodcoordinate=OCMATHET.periodcoordinate;
    out.solverinfo.monodromy=OCMATCONT.monodromy;
end
out.solverinfo.parametervaluecoord=OCMATHET.parametervaluecoord;
out.solverinfo.parameterindex=OCMATHET.parameterindex;
out.solverinfo.solutionindex=OCMATHET.solutionindex;
out.solverinfo.limitcycleindex=OCMATHET.limitcycleindex;
out.solverinfo.equilibriumcounter=OCMATHET.equilibriumcounter;
out.solverinfo.limitcyclecounter=OCMATHET.limitcyclecounter;
out.solverinfo.limitsettype=OCMATHET.limitsettype;
out.solverinfo.octrajectory2limset=OCMATHET.octrajectory2limset;



%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATHET
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
global OCMATHET OCBVP

OCBVP.limtcycleindex=OCMATHET.limitcycleindex;

% ------------------------------------------------------

function WorkspaceDone


% ------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATHET

modelpar=OCMATHET.parametervalue;
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
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATHET OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATHET=OCMATHET;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATHET.basicglobalvarfilename '4heteroclinicep2lc'],'MODELINFO')
    end
    save([OCMATHET.basicresultfilename '4heteroclinicep2lc'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATHET

pathname=OCMATHET.datapath();


%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATHET

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);


for ii=1:OCMATHET.hetorder
    if ~OCMATHET.simple(ii)
        Y=freepar(OCMATHET.Ycoordinate{ii});
        OCMATCONT.adapted = 1;
        %
        [U,S,V]=svd(OCMATHET.Q0{ii}(:,1:OCMATHET.subspacedim{ii})+OCMATHET.Q0{ii}(:,OCMATHET.subspacedim{ii}+1:end)*Y);
        OCMATHET.Q0{ii}= U;
        OCMATHET.Y{ii}=zeros(OCMATHET.orthspacedim{ii},OCMATHET.subspacedim{ii});

        freepar(OCMATHET.Ycoordinate{ii})=OCMATHET.Y{ii};
    elseif strcmp(OCMATHET.limitsettype{ii},'e')
        hatx=freepar(OCMATHET.equilibriumcoordinate{OCMATHET.octrajectory2limset(ii,3)});
        Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,OCMATHET.limsetarcarg{ii});
        if OCMATHET.implicit && OCMATHET.implicitcontrolnum(ii)
            dudx=OCMATHET.dimplicitcontroldx(0,hatx,modelpar,OCMATHET.limsetarcarg{ii});
            Jac=Jac(OCMATHET.statecostatecoordinate,OCMATHET.statecostatecoordinate)+Jac(OCMATHET.statecostatecoordinate,OCMATHET.statecostatecoordinate(end)+(1:OCMATHET.implicitcontrolnum(ii)))*dudx;
        end
        OCMATHET.asymptoticmatrix{ii}=asymptoticbc(Jac,OCMATHET.pathtype{ii},OCMATHET.maptype{ii},OCMATCONT.ZeroDeviationTolerance(ii),OCMATCONT.AsymptoticBCMethod{ii});
    end
    if isfield(OCMATHET,'fixlccoordinate')&& ~OCMATHET.fixlccoordinate(ii)
        OCMATHET.initialpoint{ii}=y(:,1);
        OCMATHET.velocityvector{ii}=ode(0,y(:,1),1,freepar,modelpar);
        OCMATHET.velocityvector{ii}=OCMATHET.velocityvector{ii}/norm(OCMATHET.velocityvector{ii});
    end
end
coeff(OCMATCONT.HE.parametercoord)=freepar;

flag = 1;


%-----------------------------------------------------------------
function out=ricatti(Y,J,Q0,dim)
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=ricatticoefficient(Q0,J,dim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
