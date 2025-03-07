function out=indifferencesolution4emf()

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

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if ~OCMATINDIF.freeendtime(solutionindex)
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime{solutionindex}];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.truncationtimecoord{solutionindex}).'];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if ~OCMATINDIF.freeendtime(solutionindex)
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime{solutionindex}];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.truncationtimecoord{solutionindex}).'];
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
        if OCMATINDIF.freeendtime(solutionindex)
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.truncationtimecoord{solutionindex})=dxdt;
        end
    end
else
    if OCMATINDIF.freeendtime(solutionindex)
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATINDIF.truncationtimecoord{solutionindex})=dxdt;
    end
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT
resinit=[];
resasym=[];
resconnec=[];
residpt=[];
resemf=[];
resricatti=[];

for order=1:OCMATINDIF.indifferenceorder
    switchtimes=freepar(OCMATINDIF.switchtimecoord{order});
    arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{order});
    % calculate equilibrium
    hatx(OCMATINDIF.emfcoord{order},1)=freepar(OCMATINDIF.emfindex{order});
    if ~isempty(OCMATINDIF.explicitemfcoordinate{order})
        hatx(OCMATINDIF.explicitemfcoordinate{order},1)=OCMATINDIF.explicitequilibriumvalue(freepar(OCMATINDIF.emfindex{order}),modelpar,arcarg(end));
    end
    Y=freepar(OCMATINDIF.Yindex{order});
    asymptoticmatrix=OCMATINDIF.Q0{order}*[-Y';OCMATINDIF.Id{order}];
    Jac=OCMATINDIF.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
    resricatti=[resricatti;ricatti(Y,Jac,order)];
    if ~isempty(OCMATINDIF.dependentemfcoordinate{order})
        tmp=OCMATINDIF.equilibrium(hatx,modelpar,arcarg(end));
        resemf=[resemf;tmp(OCMATINDIF.dependentemfcoordinate{order})];
    end

    if order==1
        if ~isempty(OCMATINDIF.freevector)
            initialstate=OCMATINDIF.startvalue;
            for jj=1:length(OCMATINDIF.freevectorcoord)
                initialstate=initialstate+freepar(OCMATINDIF.freevectorcoord(jj))*OCMATINDIF.freevector(:,jj);
            end
            restarget=depvara(OCMATINDIF.statecoordinate,1)-initialstate;
        end
    end
    if order<OCMATINDIF.indifferenceorder
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(order+1))-depvara(OCMATINDIF.statecoordinate,OCMATINDIF.initcoord(order))];
        residpt=[residpt; ...
            OCMATINDIF.bcindifference(depvara,modelpar,OCMATCONT.HE.arcarg([OCMATINDIF.arccoord{:}]),OCMATINDIF.initcoord([order order+1]))];
    end

    resasym=[resasym;OCMATINDIF.bcasymptotic(depvarb(:,OCMATINDIF.cumsumnumarc(order)),asymptoticmatrix,hatx)];
    for arc=1:OCMATINDIF.numarc(order)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset(depvara(:,OCMATINDIF.arccoord{order}),depvarb(:,OCMATINDIF.arccoord{order}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{order},arc); ...
            OCMATINDIF.guard(depvara(:,OCMATINDIF.arccoord{order}),depvarb(:,OCMATINDIF.arccoord{order}),modelpar,switchtimes,arcarg,OCMATINDIF.edge{order},arc)];
    end
end
res=[restarget;resinit;resasym;resconnec;residpt;resemf;resricatti];

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
    if ~OCMATINDIF.freeendtime(solutionindex)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime{solutionindex}];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.truncationtimecoord{solutionindex}).'];
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
    if OCMATINDIF.freeendtime(order)==0
        hatx=OCMATINDIF.saddlepoint{order};
        yend=y(:,rightarcindex(OCMATINDIF.cumsumnumarc(order)));
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
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
for arc=1:sum(OCMATINDIF.cumsumnumarc(end))
    solutionindex=OCMATINDIF.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
    end
    if ~OCMATINDIF.freeendtime(solutionindex)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.truncationtime{solutionindex}];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.truncationtimecoord{solutionindex}).'];
    end
    diffarctime=diff(arctime);
    if solutionindex>1
        transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        transformedtimeshift=0;
    end
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(relarc)*(s(leftarcindex(arc):rightarcindex(arc))-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
end
% clear possible persistent variable
h=OCMATINDIF.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar]=drearr(tmesh,coeff,tangent);
if ~isempty(OCMATINDIF.targetvalue)
    if ~isempty(OCMATINDIF.targetcoordinate)
        fprintf(1,' Difference to target value: %g\n',OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1));
    elseif ~isempty(OCMATINDIF.targetvectorcoordinate)
        fprintf(1,' Difference to target value: %g\n',OCMATINDIF.targetvalue-freepar(OCMATINDIF.freevectorcoord(OCMATINDIF.targetvectorcoordinate)));
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
global OCMATINDIF

failed=[];
for ii=id
    switch ii
        case 1
            if ~OCMATINDIF.stopcriterion
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                if ~isempty(OCMATINDIF.targetcoordinate)
                    out=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
                elseif ~isempty(OCMATINDIF.targetvectorcoordinate)
                    out=OCMATINDIF.targetvalue-freepar(OCMATINDIF.freevectorcoord(OCMATINDIF.targetvectorcoordinate));
                end
            else
                out=tangent(end);
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF OCBVP
[t,y,z,freepar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,OCMATINDIF.parametervalue);
out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
    OCMATCONT.HE.TIMEDDATA.rightarcindex];
out.arcinterval=[];
out.timehorizon=[];
for ii=1:OCMATINDIF.indifferenceorder
    arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii});

    if ~OCMATINDIF.freeendtime(ii)
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' OCMATINDIF.truncationtime{ii}];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{ii}).' freepar(OCMATINDIF.truncationtimecoord{ii}).'];
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
    hatx=OCMATINDIF.saddlepoint{ii};
    hatx(OCMATINDIF.emfcoord{ii},1)=freepar(OCMATINDIF.emfindex{ii});
    if ~isempty(OCMATINDIF.explicitemfcoordinate{ii})
        hatx(OCMATINDIF.explicitemfcoordinate{ii},1)=OCMATINDIF.explicitequilibriumvalue(freepar(OCMATINDIF.emfindex{ii}),modelpar,arcarg(end));
    end
    out.solverinfo.saddlepoint{ii}=hatx;
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATINDIF.initialtime;
out.modelparameter=OCMATINDIF.parametervalue;
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
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.truncationtimecoord=OCMATINDIF.truncationtimecoord;
out.solverinfo.emfcoord=OCMATINDIF.emfcoord;
out.solverinfo.emfindex=OCMATINDIF.emfindex;
out.solverinfo.Yindex=OCMATINDIF.Yindex;
out.solverinfo.emfindex=OCMATINDIF.emfindex;
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

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATINDIF.parametervalue;
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
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
    otherwise
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
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4emf'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4emf'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATINDIF

pathname=OCMATINDIF.datapath();

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATAE

discretizationdata=OCMATAE.domaindiscretization(arcarg);

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
function out=ricatti(Y,J,idx)
global OCMATINDIF
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(OCMATINDIF.Q0{idx},J,OCMATINDIF.subspacedim(idx));
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
