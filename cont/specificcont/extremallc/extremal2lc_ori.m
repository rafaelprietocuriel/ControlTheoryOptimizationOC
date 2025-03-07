function out=extremal2lc()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;
out{5}{3}=@odehess;
out{5}{4}=@bchess;
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
global OCMATCONT

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
if OCMATCONT.monodromy
    [J,M]=calc_RHSJac4LC(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
    OCMATCONT.monodromymatrix=M;
else
    J=calc_RHSJac(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
end
% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,icfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end
% J=Jnum;

function opt=options(opt)

function varargout=probleminit(varargin)
tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;

function [F,J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATCONT OCMATAE

switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.continuationparameterindex)=freemodelparameter;
end
solutionindex=OCMATAE.solutionindex(arc);
if solutionindex==1
    relarc=arc;
    transformedtimeshift=0;
    if strcmp(OCMATCONT.continuationtype,'time')
        arctime=[OCMATAE.TRJ.switchtimes OCMATAE.TRJ.endtime];
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.TRJ.freetimecoordinate);
        arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        arctime=[0 arctime];
        if OCMATAE.freeendtime
            arctime(end)=freepar(OCMATAE.TRJ.endtimecoordinate);
        end
    else
        if OCMATAE.freeendtime
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' freepar(OCMATAE.TRJ.endtimecoordinate).'];
        else
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' OCMATAE.TRJ.endtime];
        end
    end
    arcarg=OCMATAE.TRJ.arcarg(arc);
else
    relarc=arc-OCMATAE.TRJ.numarc;
    transformedtimeshift=OCMATAE.TRJ.numarc;
    arctime=[0 freepar(OCMATAE.LC.switchtimecoordinate).' freepar(OCMATAE.LC.periodcoordinate).'];
    arcarg=OCMATAE.LC.arcarg(relarc);
end

diffarctime=diff(arctime);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
solutionindex=OCMATAE.solutionindex(arc);
if solutionindex==1
    relarc=arc;
    transformedtimeshift=0;
    if strcmp(OCMATCONT.continuationtype,'time')
        arctime=[OCMATAE.TRJ.switchtimes OCMATAE.TRJ.endtime];
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.TRJ.freetimecoordinate);
        arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        arctime=[0 arctime];
        if OCMATAE.freeendtime
            arctime(end)=freepar(OCMATAE.TRJ.endtimecoordinate);
        end
    else
        if OCMATAE.freeendtime
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' freepar(OCMATAE.TRJ.endtimecoordinate).'];
        else
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' OCMATAE.TRJ.endtime];
        end
    end
    arcarg=OCMATAE.TRJ.arcarg(arc);
else
    relarc=arc-OCMATAE.TRJ.numarc;
    transformedtimeshift=OCMATAE.TRJ.numarc;
    arctime=[0 freepar(OCMATAE.LC.switchtimecoordinate).' freepar(OCMATAE.LC.periodcoordinate).'];
    arcarg=OCMATAE.LC.arcarg(relarc);
end

diffarctime=diff(arctime);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;

Jpar=zeros(OCBVP.numode,OCMATCONT.HE.numparameter);
dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
if solutionindex==1
    if OCMATAE.TRJ.numarc>1
        if ~OCMATAE.autonomous
            Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
        else
            Jt=0;
        end
        if relarc<OCMATAE.TRJ.numarc
            if any(relarc==OCMATAE.freetimeindex)
                Jpar(1:OCBVP.numode,OCMATAE.TRJ.switchtimecoordinate(relarc))=dxdt+diffarctime(relarc)*(s-relarc+1)*Jt;
            end
            if relarc>1
                if any(relarc-1==OCMATAE.freetimeindex)
                    Jpar(1:OCBVP.numode,OCMATAE.TRJ.switchtimecoordinate(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
                end
            end
        else
            if any(relarc==OCMATAE.freetimeindex)
                Jpar(1:OCBVP.numode,OCMATAE.TRJ.switchtimecoordinate(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
            end
            if OCMATAE.freeendtime
                Jpar(1:OCBVP.numode,OCMATAE.TRJ.endtimecoordinate)=dxdt;
                if any(relarc-1==OCMATAE.freetimeindex)
                    Jpar(1:OCBVP.numode,OCMATAE.TRJ.switchtimecoordinate(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
                end
            end
        end
    else
        if OCMATAE.freeendtime
            Jpar(1:OCBVP.numode,OCMATAE.TRJ.endtimecoordinate)=dxdt;
        end
    end
else
    if OCMATAE.LC.numarc>1
        if ~OCMATAE.autonomous
            Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
        else
            Jt=0;
        end
        if relarc<OCMATAE.LC.numarc
            Jpar(1:OCBVP.numode,OCMATAE.LC.switchtimecoordinate(relarc))=dxdt+diffarctime(relarc)*(s-relarc+1)*Jt;
            if relarc>1
                Jpar(1:OCBVP.numode,OCMATAE.LC.switchtimecoordinate(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
            end
        else
            Jpar(1:OCBVP.numode,OCMATAE.LC.switchtimecoordinate(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
            Jpar(1:OCBVP.numode,OCMATAE.LC.periodcoordinate)=dxdt;
        end
    else
        Jpar(1:OCBVP.numode,OCMATAE.LC.periodcoordinate)=dxdt;
    end
end
switch OCMATCONT.continuationtype
    case 'parameter'
        Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
        if isempty(OCMATAE.targetvalue)
            Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATAE.continuationparameterindex);
        else
            Jpar(:,OCMATCONT.HE.numparameter)=OCMATAE.continuationvector.'*Jmodelpar(:,OCMATAE.continuationparameterindex);
        end
    case 'time'
        if ~OCMATAE.autonomous
            Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
        else
            Jt=0;
        end
        if relarc==OCMATAE.continuationtimeindex
            Jpar(1:OCBVP.numode,OCMATAE.continuationcoordinate)=dxdt+diffarctime(relarc)*(s-relarc+1)*Jt;
        elseif any(relarc-1==OCMATAE.continuationtimeindex)
            Jpar(1:OCBVP.numode,OCMATAE.continuationcoordinate)=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [Hy2,Hypar,Hpar2]=odehess(s,depvar,arc,freepar,modelpar)
Hy2=[];
Hypar=[];
Hpar2=[];

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT
resconnec=[];
resricatti=[];
resper=[];

switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.continuationparameterindex)=freemodelparameter;
        switchtimes=freepar(OCMATAE.LC.switchtimecoordinate);
        arcarg=OCMATAE.LC.arcarg;

        Y=freepar(OCMATAE.Ycoordinate);
        initialstate=OCMATAE.TRJ.initialpoint(OCMATAE.statecoordinate);
        asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
        resricatti=ricatti(Y,OCMATCONT.monodromymatrix);
        initialcoordinate=OCMATAE.fixinitstate;
        resper=OCMATAE.bclimitcycle(depvara(:,OCMATAE.limitcycleindex),depvarb(:,OCMATAE.limitcycleindex));
        resper=[resper; ...
            sum(OCMATAE.LC.velocityvector.*(depvara(OCMATAE.LC.velocitycoordinate,OCMATAE.limitcycleindex(1))-OCMATAE.LC.initialpoint))];
        saddlepoint=depvara(:,OCMATAE.limitcycleindex(1));
        for arc=1:OCMATAE.LC.numarc-1
            resconnec=[resconnec; ...
                OCMATAE.reset(depvara(:,OCMATAE.limitcycleindex),depvarb(:,OCMATAE.limitcycleindex),modelpar,switchtimes,arcarg,OCMATAE.LC.edge,arc); ...
                OCMATAE.guard(depvara(:,OCMATAE.limitcycleindex),depvarb(:,OCMATAE.limitcycleindex),modelpar,switchtimes,arcarg,OCMATAE.LC.edge,arc)];
        end
        switchtimes=freepar(OCMATAE.TRJ.switchtimecoordinate);

    case 'initialstate'
        initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
        initialcoordinate=OCMATAE.TRJ.continuationcoordinate;
        asymptoticmatrix=OCMATAE.asymptoticmatrix;
        saddlepoint=OCMATAE.LC.initialpoint;
        switchtimes=freepar(OCMATAE.TRJ.switchtimecoordinate);
    case 'time'
        switchtimes=OCMATAE.TRJ.switchtimes;
        switchtimes(OCMATAE.freetimeindex)=freepar(OCMATAE.TRJ.freetimecoordinate);
        switchtimes(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        initialstate=OCMATAE.TRJ.initialpoint(OCMATAE.statecoordinate);
        initialcoordinate=OCMATAE.fixinitstate;
        asymptoticmatrix=OCMATAE.asymptoticmatrix;
        saddlepoint=OCMATAE.LC.initialpoint;
end

resinit=OCMATAE.bcinitial(depvara(:,OCMATAE.trajectoryindex),initialcoordinate,initialstate);
arcarg=OCMATAE.TRJ.arcarg;

resasym=OCMATAE.bcasymptotic(depvarb(:,OCMATAE.trajectoryindex(end)),asymptoticmatrix,saddlepoint);
if OCMATAE.fixdistance
    resasym=[resasym; ...
        sqrt(sum((saddlepoint-depvarb(:,OCMATAE.trajectoryindex(end))).^2))-OCMATAE.distance];
end
for arc=1:OCMATAE.TRJ.numarc-1
    if any(arc==OCMATAE.freetimeindex)
        resconnec=[resconnec; ...
            OCMATAE.reset(depvara(:,OCMATAE.trajectoryindex),depvarb(:,OCMATAE.trajectoryindex),modelpar,switchtimes,arcarg,OCMATAE.TRJ.edge,arc); ...
            OCMATAE.guard(depvara(:,OCMATAE.trajectoryindex),depvarb(:,OCMATAE.trajectoryindex),modelpar,switchtimes,arcarg,OCMATAE.TRJ.edge,arc)];
    else
        resconnec=[resconnec; ...
            OCMATAE.reset(depvara(:,OCMATAE.trajectoryindex),depvarb(:,OCMATAE.trajectoryindex),modelpar,switchtimes,arcarg,OCMATAE.TRJ.edge,arc)];
    end
end

res=[resinit;resconnec;resasym;resricatti;resper];

%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATAE.solutionindex(arc);
    if solutionindex==1
        relarc=arc;
        transformedtimeshift=0;
        if strcmp(OCMATCONT.continuationtype,'time')
            arctime=[OCMATAE.TRJ.switchtimes OCMATAE.TRJ.endtime];
            arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.TRJ.freetimecoordinate);
            arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
            arctime=[0 arctime];
            if OCMATAE.freeendtime
                arctime(end)=freepar(OCMATAE.TRJ.endtimecoordinate);
            end
        else
            if OCMATAE.freeendtime
                arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' freepar(OCMATAE.TRJ.endtimecoordinate).'];
            else
                arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' OCMATAE.TRJ.endtime];
            end
        end
        arcarg=OCMATAE.TRJ.arcarg(arc);
    else
        relarc=arc-OCMATAE.TRJ.numarc;
        transformedtimeshift=OCMATAE.TRJ.numarc;
        arctime=[0 freepar(OCMATAE.LC.switchtimecoordinate).' freepar(OCMATAE.LC.periodcoordinate).'];
        arcarg=OCMATAE.LC.arcarg(relarc);
    end

    diffarctime=diff(arctime);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATAE.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        [constr,labelS]=OCMATAE.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    if OCMATAE.reversetime && solutionindex==1
        violationmat=-diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    else
        violationmat=diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    end
    if any(violationmat)
        counter=counter+1;
        cols=find(violationmat);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=relarc;
        infoS(counter).rows='switchtime';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=diffarctime(relarc);
        infoS(counter).minval=min(diffarctime(:));
        b=min([b infoS(counter).minval]);
    end

    switch OCMATCONT.continuationtype
        case 'parameter'
            depvara=y(:,leftarcindex);
            depvarb=y(:,rightarcindex);
            hatx=depvara(:,OCMATAE.limitcycleindex(1));
            yend=depvarb(:,OCMATAE.trajectoryindex(end));
        case 'initialstate'
            hatx=OCMATAE.LC.initialpoint;
            yend=y(:,end);
    end

    if ~OCMATAE.fixdistance
        violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx)<0;
        if violationmat
            counter=counter+1;
            cols=size(y,2);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='maxdistance';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=norm(norm(yend-hatx));
            infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx);
            b=min([b infoS(counter).minval]);
        end
    end

end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);

leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex(:,OCMATAE.solutionindex==1);
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex(:,OCMATAE.solutionindex==1);
s=s(leftarcindex(1):rightarcindex(end));
t=s;
y=y(:,leftarcindex(1):rightarcindex(end));
if strcmp(OCMATCONT.continuationtype,'time')
    arctime=[OCMATAE.TRJ.switchtimes OCMATAE.TRJ.endtime];
    arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.TRJ.freetimecoordinate);
    arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
    arctime=[0 arctime];
    if OCMATAE.freeendtime
        arctime(end)=freepar(OCMATAE.TRJ.endtimecoordinate);
    end
else
    if OCMATAE.freeendtime
        arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' freepar(OCMATAE.TRJ.endtimecoordinate).'];
    else
        arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' OCMATAE.TRJ.endtime];
    end
end
diffarctime=diff(arctime);
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATAE.solutionindex(arc);
    if solutionindex==1
        t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*(s(leftarcindex(arc):rightarcindex(arc)))+(arctime(arc)-diffarctime(arc)*(arc-1));
    end
end

h=OCMATAE.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg(OCMATAE.solutionindex==1),freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT
idx=[];
if isempty(coeff)
    return
end
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATAE
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
global OCMATCONT OCMATAE
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
global OCMATCONT

failed=[];
for ii=id
    switch ii
        case 1
            switch continuationtype
                case {'initialstate','parameter'}
                    out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
                case 'time'
                    if ~isempty(OCMATAE.targetvalue)
                        out=OCMATAE.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
                    else
                        out=[];
                    end
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATAE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end

switch OCMATCONT.continuationtype
    case 'initialstate'
        if OCMATAE.freeendtime
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' freepar(OCMATAE.TRJ.endtimecoordinate).'];
        else
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' OCMATAE.TRJ.endtime];
        end
        timehorizon=arctime(end);
    case 'parameter'
        if OCMATAE.freeendtime
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' freepar(OCMATAE.TRJ.endtimecoordinate).' 0 freepar(OCMATAE.LC.switchtimecoordinate).' freepar(OCMATAE.LC.periodcoordinate).'];
            timehorizon=[freepar(OCMATAE.TRJ.endtimecoordinate) freepar(OCMATAE.LC.periodcoordinate)];
        else
            arctime=[0 freepar(OCMATAE.TRJ.switchtimecoordinate).' OCMATAE.TRJ.endtime 0 freepar(OCMATAE.LC.switchtimecoordinate).' freepar(OCMATAE.LC.periodcoordinate).'];
            timehorizon=[OCMATAE.TRJ.endtime freepar(OCMATAE.LC.periodcoordinate)];
        end

        out.solverinfo.solutionindex=OCMATAE.solutionindex;
        out.solverinfo.Ycoordinate=OCMATAE.Ycoordinate;
        out.solverinfo.lcswitchtimecoordinate=OCMATAE.LC.switchtimecoordinate;
        out.solverinfo.lcperiodcoordinate=OCMATAE.LC.periodcoordinate;
        out.solverinfo.lcarcarg=OCMATAE.LC.arcarg;
        out.solverinfo.lcnumarc=OCMATAE.LC.numarc;
        out.solverinfo.monodromy=OCMATCONT.monodromy;
    case 'time'
        arctime=[OCMATAE.TRJ.switchtimes OCMATAE.TRJ.endtime];
        arctime(OCMATAE.freetimeindex)=freepar(OCMATAE.TRJ.freetimecoordinate);
        arctime(OCMATAE.continuationtimeindex)=freepar(OCMATAE.continuationcoordinate);
        arctime=[0 arctime];
        if OCMATAE.freeendtime
            arctime(end)=freepar(OCMATAE.TRJ.endtimecoordinate);
        end
        timehorizon=arctime(end);
end
out.arcarg=OCMATAE.TRJ.arcarg;
out.arcinterval=arctime;
out.timehorizon=timehorizon;

out.x0=0;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;

out.solverinfo.trjswitchtimecoordinate=OCMATAE.TRJ.switchtimecoordinate;
if OCMATAE.freeendtime
    out.solverinfo.trjendtimecoordinate=OCMATAE.TRJ.endtimecoordinate;
end
out.solverinfo.trjarcarg=OCMATAE.TRJ.arcarg;
out.solverinfo.trjnumarc=OCMATAE.TRJ.numarc;
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;

out.solverinfo.continuationclass='extremal2lc';
out.solverinfo.continuationtype=OCMATCONT.continuationtype;
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE

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
global OCMATAE OCBVP

OCBVP.limtcycleindex=OCMATAE.solutionindex==2;
% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATAE

z=[];
modelpar=OCMATAE.parametervalue;
y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
switch OCMATCONT.continuationtype
    case 'parameter'
        if ~isempty(OCMATAE.targetvalue)
            freemodelparameter=OCMATAE.startvalue+freepar(OCMATAE.continuationcoordinate)*OCMATAE.continuationvector;
        else
            freemodelparameter=freepar(OCMATAE.continuationcoordinate);
        end
        modelpar(OCMATAE.continuationparameterindex)=freemodelparameter;
end

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATAE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATAE=OCMATAE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATAE.basicglobalvarfilename '4extremal2lc'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremal2lc'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATAE

pathname=OCMATAE.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATAE

switch OCMATCONT.continuationtype
    case 'parameter'

        leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
        [t,y,z,freepar]=drearr(tmesh,coeff);
        Y=freepar(OCMATAE.Ycoordinate);
        OCMATCONT.adapted = 1;
        %
        [U,S,V]=svd(OCMATAE.Q0(:,1:OCMATAE.subspacedim)+OCMATAE.Q0(:,OCMATAE.subspacedim+1:end)*Y);
        OCMATAE.Q0=U;
        OCMATAE.Y=zeros(OCMATAE.orthspacedim,OCMATAE.subspacedim);

        freepar(OCMATAE.Ycoordinate)=OCMATAE.Y;
        
        depvara=y(:,leftarcindex);
        OCMATAE.initialpointlc=depvara(:,OCMATAE.limitcycleindex(1));

        switch OCMATCONT.bvpmethod
            case {'bvp6c','bvp4c'}
                coeff=[y(:);freepar];
            otherwise
        end
end
flag = 1;

%-----------------------------------------------------------------
function out=ricatti(Y,J)
global OCMATAE
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(OCMATAE.Q0,J,OCMATAE.subspacedim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
