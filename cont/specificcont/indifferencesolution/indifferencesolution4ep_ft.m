function out=indifferencesolution4ep_ft()

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
out{12}{1}=@residual;
out{12}{2}=@maxresidual;
out{12}{3}=@verifyresidual;
out{12}{4}=@prepare4meshadapt;
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
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);%OCMATINDIF.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if OCMATINDIF.movinghorizon(solutionindex)
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.movinghorizoncoord(solutionindex))];
elseif OCMATINDIF.optimalhorizon(solutionindex)
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.optimalhorizoncoord{solutionindex})];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if OCMATINDIF.movinghorizon(solutionindex)
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.movinghorizoncoord(solutionindex))];
elseif OCMATINDIF.optimalhorizon(solutionindex)
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.optimalhorizoncoord{solutionindex})];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATINDIF.objectivevaluecalc
    J=[J; ...
        dtds*OCMATINDIF.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
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
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc))=dxdt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-dxdt;
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        if OCMATINDIF.movinghorizon(solutionindex)
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.truncationtimecoord{solutionindex})=dxdt;
        elseif OCMATINDIF.optimalhorizon(solutionindex)
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.optimalhorizoncoord{solutionindex})=dxdt;
        end
    end
else
    if OCMATINDIF.movinghorizon(solutionindex)
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATINDIF.truncationtimecoord{solutionindex})=dxdt;
    elseif OCMATINDIF.optimalhorizon(solutionindex)
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATINDIF.optimalhorizoncoord{solutionindex})=dxdt;
    end
end
Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATINDIF.parametervaluecoord)=Jmodelpar(:,OCMATINDIF.parameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP

resinit=[];
resasym=[];
restrans=[];
resconnec=[];
resequilibrium=[];
resricatti=[];
restarget=[];
residpt=[];

modelpar(OCMATINDIF.parameterindex)=freepar(OCMATINDIF.parametervaluecoord);

for ii=1:OCMATINDIF.indifferenceorder
    if ii==1
        if ~isempty(OCMATINDIF.freevector)
            initialstate=OCMATINDIF.startvalue;
            for jj=1:length(OCMATINDIF.freevectorcoord)
                initialstate=initialstate+freepar(OCMATINDIF.freevectorcoord(jj))*OCMATINDIF.freevector(:,jj);
            end
            restarget=depvara(OCMATINDIF.statecoordinate,1)-initialstate;
        end
    end
    switchtimes=freepar(OCMATINDIF.switchtimecoord{ii});
    arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii});
    if ii<OCMATINDIF.indifferenceorder
        if any(strcmp(OCMATINDIF.pathtype(ii:ii+1),'sts'))
            resinit=[resinit; ...
                depvara(OCMATINDIF.statecoordinate(1),OCMATINDIF.initcoord(ii+1))-depvara(OCMATINDIF.statecoordinate(1),OCMATINDIF.initcoord(ii))];
        else
            resinit=[resinit; ...
                depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(ii+1))-depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(ii))];
        end
        residpt=[residpt; ...
            OCMATINDIF.bcindifference(depvara,modelpar,OCMATCONT.HE.arcarg([OCMATINDIF.arccoord{:}]),OCMATINDIF.initcoord([ii ii+1]))];
    end
    if OCMATINDIF.ocasymptotic(ii)
        hatx=freepar(OCMATCONT.HE.equilibriumcoord{ii});
        resequilibrium=[resequilibrium;OCMATINDIF.canonicalsystem(0,hatx,modelpar,arcarg(end))];
        if OCMATINDIF.simple
            asymptoticmatrix=OCMATINDIF.asymptoticmatrix{ii};
        else
            Y=freepar(OCMATCONT.HE.Ycoord{ii});
            asymptoticmatrix=OCMATINDIF.Q0{ii}*[-Y';OCMATINDIF.Id{ii}];
            Jac=OCMATINDIF.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
            resricatti=[resricatti;ricatti(Y,Jac,OCMATINDIF.Q0{ii},OCMATINDIF.subspacedim{ii})];
        end
        resasym=[resasym;OCMATINDIF.bcasymptotic(depvarb(:,OCMATINDIF.cumsumnumarc(ii)),asymptoticmatrix,hatx)];
        for arc=1:OCMATINDIF.numarc(ii)-1
            resconnec=[resconnec; ...
                OCMATINDIF.reset(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
                OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
        end
        if OCMATINDIF.movinghorizon(ii)
            resasym=[resasym; ...
                sqrt(sum((hatx-depvarb(:,OCMATINDIF.cumsumnumarc(ii))).^2))-OCMATINDIF.distance(ii)];
        end
    else
        if ~OCMATINDIF.optimalhorizon(ii)
            timehorizon=OCMATINDIF.truncationtime(ii);
        else
            timehorizon=freepar(OCMATINDIF.optimalhorizoncoord{ii});
        end
        restrans=OCMATINDIF.bctransversality(timehorizon,depvarb(:,OCMATINDIF.cumsumnumarc(ii)),modelpar,arcarg(end));
        if ~isempty(OCMATINDIF.fixendstatecoord{ii})
            restrans(OCMATINDIF.fixendstatecoord{ii},1)=depvarb(OCMATINDIF.fixendstatecoord{ii},OCMATINDIF.cumsumnumarc(ii))-OCMATINDIF.endstate{ii};
        end
        if ~isempty(OCMATINDIF.fixinitstatecoord{ii})
            resinit=[resinit; ...
                depvara(OCMATINDIF.fixinitstatecoord{ii},OCMATINDIF.initcoord(ii))-OCMATINDIF.initstate{ii}];
        end
        if OCMATINDIF.optimalhorizon(ii)
            restrans=[restrans; ...
                OCMATINDIF.bcoptimalhorizon(timehorizon,depvarb(:,OCMATINDIF.cumsumnumarc(ii)),modelpar,arcarg(end))];
        end
        for arc=1:OCMATINDIF.numarc(ii)-1
            resconnec=[resconnec; ...
                OCMATINDIF.reset(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
                OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
        end
    end
end
res=[resinit;resasym;resconnec;resequilibrium;restrans;resricatti;restarget;residpt];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
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
    solutionindex=OCMATINDIF.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    if OCMATINDIF.movinghorizon(solutionindex)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.movinghorizoncoord(solutionindex))];
    elseif OCMATINDIF.optimalhorizon(solutionindex)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.optimalhorizoncoord{solutionindex})];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime(solutionindex)];
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
        [constr labelS]=OCMATINDIF.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
%     if OCMATINDIF.stableflag{arc}
%         violationmat=diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
%     else
%         violationmat=-diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
%     end
%     if any(violationmat)
%         counter=counter+1;
%         cols=find(violationmat);
%         infoS(counter).arcarg=arcarg;
%         infoS(counter).arcnum=arc;
%         infoS(counter).rows='switchtime';
%         infoS(counter).cols=cols;
%         infoS(counter).violationmat=violationmat;
%         %infoS(counter).constraintvalue=diffarctime(arc);
%         infoS(counter).minval=min(diffarctime(:));
%         b=min([b infoS(counter).minval]);
%     end
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
h=OCMATINDIF.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
if ~isempty(OCMATINDIF.hitstatevalue)
    fprintf(1,' Difference hitvalue: %g\n',y(OCMATINDIF.hitstatecoordinate,1)-OCMATINDIF.hitstatevalue);
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

failed=[];
for ii=id
    switch ii
        case 1
            if ~isempty(OCMATINDIF.hitstatevalue)
                out=y(OCMATINDIF.hitstatecoordinate,1)-OCMATINDIF.hitstatevalue;
            else
                out=OCMATINDIF.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
            end
            %out=tangent(end);
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
timehorizon=zeros(1,OCMATINDIF.indifferenceorder);
for ii=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.movinghorizon(ii)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' freepar(OCMATINDIF.movinghorizoncoord(ii))];
    elseif OCMATINDIF.optimalhorizon(ii)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' freepar(OCMATINDIF.optimalhorizoncoord{ii})];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' OCMATINDIF.truncationtime(ii)];
    end
    out.arcinterval=[out.arcinterval arctime];
    timehorizon(ii)=arctime(end);
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATINDIF.initialtime;
out.timehorizon=timehorizon;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution4ep_ft';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATINDIF.inftimetransformation;
out.solverinfo.pathtype=OCMATINDIF.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.truncationtimecoord=OCMATINDIF.truncationtimecoord;
out.solverinfo.equilibriumcoord=OCMATCONT.HE.equilibriumcoord;
out.solverinfo.Ycoord=OCMATCONT.HE.Ycoord;
out.solverinfo.subspacedim=OCMATINDIF.subspacedim;
out.solverinfo.orthspacedim=OCMATINDIF.orthspacedim;
out.solverinfo.qbasis=OCMATINDIF.Q0;
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

%numarc=OCMATCONT.HE.numarc;
%domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATINDIF.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
modelpar(OCMATINDIF.parameterindex)=freepar(OCMATINDIF.parametervaluecoord);

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATINDIF OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATINDIF=OCMATINDIF;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4ep_ft'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4ep_ft'],'sout','bvpout')
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
        warning('off','MATLAB:deval:NonuniqueSolution');
        ynew=devalbvpoc(sol,tmeshnew);
        warning('on','MATLAB:deval:NonuniqueSolution');
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,n,N,nN);

        coeffnew=[ynew(:);freepar(:)];
end



%-----------------------------------------------------------------
function out=ricatti(Y,J,Q0,dim)
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(Q0,J,dim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
