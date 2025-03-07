    function out=indifferencesolutionp()

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
modelpar(OCMATINDIF.parameterindex)=freepar(OCMATINDIF.parametervaluecoordinate);
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
%arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.numarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);

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
Jpar=zeros(numode,OCMATCONT.HE.numparameter);
if OCMATINDIF.numarc(solutionindex)>1
    dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATINDIF.autonomous
        Jt=OCMATINDIF.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(1:numode,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc))=dxdt;
        if relarc>1
            Jpar(1:numode,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc-1))=-dxdt;
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        if OCMATINDIF.freeendtime(solutionindex)
            %dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.truncationtimecoordinate(solutionindex))=dxdt;
        end
    end
else
    if OCMATINDIF.freeendtime(solutionindex)
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATINDIF.truncationtimecoordinate(solutionindex))=dxdt;
    end
end
Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATINDIF.parametervaluecoordinate)=Jmodelpar(:,OCMATINDIF.parameterindex);
if ~isempty(OCMATINDIF.freeparameter)
    Jpar(:,OCMATINDIF.freeparametercoordinate)=Jmodelpar(:,OCMATINDIF.freeparameter);
end


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT
resinit=[];
resasym=[];
resconnec=[];
residpt=[];
resequilibrium=[];
resricatti=[];
restarget=[];
resper=[];

modelpar(OCMATINDIF.parameterindex)=freepar(OCMATINDIF.parametervaluecoordinate);
if ~isempty(OCMATINDIF.freeparameter)
    modelpar(OCMATINDIF.freeparameter)=freepar(OCMATINDIF.freeparametercoordinate);
end

if OCMATINDIF.limitcyclecounter
    % conditions for limitcycle
    initaldepvarlc=zeros(OCMATINDIF.statecostatecoordinate(end),OCMATINDIF.limitcyclecounter);
    for ii=1:OCMATINDIF.limitcyclecounter
        actdepvara=depvara(:,OCMATINDIF.limitcycleindex==ii);
        actdepvarb=depvarb(:,OCMATINDIF.limitcycleindex==ii);
        switchtimes=freepar(OCMATINDIF.switchtimecoordinate{OCMATINDIF.indifforder+ii});
        arcarg=OCMATINDIF.arcarg{OCMATINDIF.indifforder+ii};
        initaldepvarlc(:,ii)=actdepvara(:,1);
        resper=[resper; ...
            OCMATINDIF.bclimitcycle(actdepvara,actdepvarb); ...
            sum(OCMATINDIF.velocityvector{ii}.*(actdepvara(OCMATINDIF.velocitycoordinate{ii},1)-OCMATINDIF.initialpoint{ii}))];

        for arc=1:OCMATINDIF.numarc(OCMATINDIF.indifforder+ii)-1
            resconnec=[resconnec; ...
                OCMATINDIF.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{OCMATINDIF.indifforder+ii},arc); ...
                OCMATINDIF.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{OCMATINDIF.indifforder+ii},arc)];
        end

    end
end
initialdepvar=zeros(size(depvara,1),OCMATINDIF.indifferenceorder);
initialarcarg=zeros(1,OCMATINDIF.indifferenceorder);
userbcvalue=[];
for ii=1:OCMATINDIF.indifferenceorder
    solutionindex=find(OCMATINDIF.solutionindex==ii);
    actdepvara=depvara(:,solutionindex);
    actdepvarb=depvarb(:,solutionindex);
    arcarg=OCMATINDIF.arcarg{ii};
    switchtimes=freepar(OCMATINDIF.switchtimecoordinate{ii});
    initialdepvar(:,ii)=actdepvara(:,1);
    initialarcarg(ii)=arcarg(1);
    if ii==1
        initialstate=OCMATINDIF.startvalue;
        if ~isempty(OCMATINDIF.freevectorcoordinate)
            for jj=1:length(OCMATINDIF.freevectorcoordinate)
                initialstate=initialstate+freepar(OCMATINDIF.freevectorcoordinate(jj))*OCMATINDIF.freevector(:,jj);
            end
            targetcoord=setdiff(OCMATINDIF.statecoordinate,OCMATINDIF.fixinitstate);
            restarget=depvara(targetcoord,1)-initialstate(targetcoord);
        end
    end
    if ii==1 & ~isempty(OCMATINDIF.fixinitstate)
        resinit=OCMATINDIF.bcinitial(initialdepvar(:,1),OCMATINDIF.fixinitstate,OCMATINDIF.startvalue(OCMATINDIF.fixinitstate),modelpar,arcarg(1));
    end
    if ii>1
        resinit=[resinit; ...
            OCMATINDIF.bcinitial(initialdepvar(:,ii),OCMATINDIF.statecoordinate,initialdepvar(OCMATINDIF.statecoordinate,ii-1),modelpar,arcarg(1))];
        residpt=[residpt; ...
            OCMATINDIF.bcindifference(initialdepvar,modelpar,initialarcarg,[ii-1 ii])];
    end
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
    if OCMATINDIF.userbc
        userbcvalue=[userbcvalue; ...
            OCMATINDIF.userbcfunc(actdepvara,actdepvarb,modelpar,arcarg)];
        if ii==2
            residpt=[residpt; ...
                userbcvalue-(1-freepar(end))*OCMATINDIF.userbcvalue];
        end
    end
    if OCMATINDIF.equilibriumcounter>1 || ii==1
        resequilibrium=[resequilibrium;OCMATINDIF.equilibrium(hatx,modelpar,arcarg(end))];
    end
    resasym=[resasym;OCMATINDIF.bcasymptotic(actdepvarb(:,end),asymptoticmatrix,hatx)];
    if OCMATINDIF.freeendtime(ii)>0
        resasym=[resasym; ...
            sqrt(sum((OCMATINDIF.saddlepoint{ii}-actdepvarb(:,end)).^2))-OCMATINDIF.distance(ii)];
    end
    for arc=1:OCMATINDIF.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
            OCMATINDIF.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
    end
end
res=[resinit;restarget;resasym;resconnec;residpt;resequilibrium;resricatti;resper];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
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
for order=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.fixdistance(order)==0
        hatx=freepar(OCMATINDIF.equilibriumcoordinate{order});
        hatx=hatx(OCMATINDIF.statecostatecoordinate);
        yend=y(OCMATINDIF.statecostatecoordinate,rightarcindex(OCMATINDIF.cumsumnumarc(order)));
        violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx)<0;
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
%sol=evalatmesh(tmesh,y,z);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
ctr=0;
clf
t=[];
for ii=1:OCMATINDIF.indifforder
    solutionindex=find(OCMATINDIF.solutionindex==ii);
    limitcycleindex=OCMATINDIF.limitcycleindex(solutionindex(1));
    if ~limitcycleindex
        ctr=ctr+1;
        transformedtimeshift=OCMATINDIF.cumsumnumarc(ctr);
        if OCMATINDIF.freeendtime(ii)
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ii}).' freepar(OCMATINDIF.truncationtimecoordinate(ii))];
        else
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ii}).' OCMATINDIF.truncationtime(ii)];
        end
        diffarctime=diff(arctime);
        for arc=1:OCMATINDIF.numarc(ii)
            t=[t diffarctime(arc)*(s(leftarcindex(solutionindex(arc)):rightarcindex(solutionindex(arc)))-transformedtimeshift)+(arctime(arc)-diffarctime(arc)*(arc-1))];
        end
    end
end
h=OCMATINDIF.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar]=drearr(tmesh,coeff,tangent);
if ~isempty(OCMATINDIF.targetparametervalue)
    out=OCMATINDIF.targetparametervalue-freepar(end);
    txt='Difference to target parameter value';
elseif ~isempty(OCMATINDIF.targetcoordinate)
    out=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
    txt='Difference to target value';
else
    out=1-freepar(end);
    txt='Continuation parameter';
end
fprintf(1,' %s: %g\n',txt,out);

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
for ii=id
    switch ii
        case 1
            [t,y,z,freepar]=drearr(tmesh,coeff,tangent);
            if ~isempty(OCMATINDIF.targetparametervalue)
                out=OCMATINDIF.targetparametervalue-freepar(end);
            elseif ~isempty(OCMATINDIF.targetcoordinate)
                out=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
            else
                out=1-freepar(end);
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

out.x=t;
out.y=y;
out.parameters=freepar;
%out=transform2nativematlab(t,coeff,OCMATINDIF.parametervalue);
out.arcinterval=[];
out.timehorizon=[];
out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
    OCMATCONT.HE.TIMEDDATA.rightarcindex];

for ii=1:OCMATINDIF.indifforder
    if OCMATINDIF.freeendtime(ii)
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ii}).' freepar(OCMATINDIF.truncationtimecoordinate(ii))];
    else
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{ii}).' OCMATINDIF.truncationtime(ii)];
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
    if ~OCMATINDIF.simple(ii)
        out.solverinfo.Ycoordinate=OCMATINDIF.Ycoordinate;
        out.solverinfo.subspacedim=OCMATINDIF.subspacedim;
        out.solverinfo.orthspacedim=OCMATINDIF.orthspacedim;
        out.solverinfo.qbasis=OCMATINDIF.Q0;
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
out.solverinfo.arcarg=OCMATINDIF.arcarg;
out.solverinfo.conttype='indifferencesolutionp';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATINDIF.pathtype;
out.solverinfo.switchtimecoordinate=OCMATINDIF.switchtimecoordinate;

out.solverinfo.equilibriumcoordinate=OCMATINDIF.equilibriumcoordinate;
out.solverinfo.Ycoordinate=OCMATINDIF.Ycoordinate;
out.solverinfo.subspacedim=OCMATINDIF.subspacedim;
out.solverinfo.orthspacedim=OCMATINDIF.orthspacedim;
out.solverinfo.qbasis=OCMATINDIF.Q0;
out.solverinfo.freevector=OCMATINDIF.freevector;
out.solverinfo.freevectorcoordinate=OCMATINDIF.freevectorcoordinate;
out.solverinfo.solutionindex=OCMATINDIF.solutionindex;
out.solverinfo.limitsettype=OCMATINDIF.limitsettype;
out.solverinfo.octrajectory2limset=OCMATINDIF.octrajectory2limset;

out.solverinfo.limitcycleindex=OCMATINDIF.limitcycleindex;
out.solverinfo.limitcyclecounter=OCMATINDIF.limitcyclecounter;
if OCMATINDIF.limitcyclecounter
    out.solverinfo.periodcoordinate=OCMATINDIF.periodcoordinate;
end
if OCMATINDIF.freeparameter
    out.solverinfo.freeparameter=OCMATINDIF.freeparameter;
    out.solverinfo.freeparametercoordinate=OCMATINDIF.freeparametercoordinate;
end
out.solverinfo.userbc=OCMATINDIF.userbc;


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
modelpar(OCMATINDIF.parameterindex)=freepar(OCMATINDIF.parametervaluecoordinate);
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
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolutionp'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolutionp'],'sout','bvpout')
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


for ii=1:OCMATINDIF.indifferenceorder
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    Y=freepar(OCMATINDIF.Ycoordinate{ii});
    OCMATCONT.adapted = 1;
    %
    [U,S,V]=svd(OCMATINDIF.Q0{ii}(:,1:OCMATINDIF.subspacedim{ii})+OCMATINDIF.Q0{ii}(:,OCMATINDIF.subspacedim{ii}+1:end)*Y);
    OCMATINDIF.Q0{ii}= U;
    OCMATINDIF.Y{ii}=zeros(OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});

    freepar(OCMATINDIF.Ycoordinate{ii})=OCMATINDIF.Y{ii};
    switch OCMATCONT.bvpmethod
        case {'bvp6c','bvp4c'}
            coeff=[y(:);freepar];
        otherwise
    end
end
flag = 1;

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

%-----------------------------------------------------------------
function out=ricatti(Y,J,Q0,dim)
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=ricatticoefficient(Q0,J,dim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
