function out=extremalpF4ft()

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
modelpar(OCMATFTE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);

depvarfund=depvar(OCMATFTE.fundcoord,:);
depvar=depvar(OCMATFTE.solcoord,:);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
dxdtfund=zeros(OCMATFTE.fundnumcoord,length(s));
for ii=1:length(s)
    depvartmp=depvarfund(:,ii);
    tmp=dtds*OCMATFTE.canonicalsystemjacobian(t,depvartmp(OCMATFTE.fundmatcoord),modelpar,arcarg);
    dxdtfund(:,ii)=tmp(:);
end
dxdt=[dxdt;dxdtfund];
if OCMATFTE.objectivevaluecalc
    dxdt(OCMATFTE.objectivevaluecoord,:)=dtds*OCMATFTE.objectivefunction(t,depvar,modelpar,arcarg);
end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
J=[];
Jpar=[];
return
modelpar(OCMATFTE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATFTE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATFTE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
Jmodelpar=dtds*OCMATFTE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
if OCMATCONT.HE.numarc>1
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
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATFTE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATFTE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATFTE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        if OCMATFTE.optimalhorizon
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATFTE.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc)*Jt;
        end
    end
else
    if ~OCMATFTE.autonomous
        Jt=OCMATFTE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATFTE.optimalhorizon
        dxdt=OCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATFTE.optimalhorizoncoord)=dxdt+diffarctime(arc)*(s-arc)*Jt;
    end
end
Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATFTE.varyparameterindex);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE OCMATCONT OCBVP
modelpar(OCMATFTE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
switchtimes=freepar(OCMATFTE.switchtimecoord);
if ~OCMATFTE.optimalhorizon
    timehorizon=OCMATFTE.truncationtime;
else
    timehorizon=freepar(OCMATFTE.optimalhorizoncoord);
end
resconnec=[];

depvarafund=depvara(OCMATFTE.fundcoord,:);
depvarbfund=depvarb(OCMATFTE.fundcoord,:);
depvara=depvara(OCMATFTE.solcoord,:);
depvarb=depvarb(OCMATFTE.solcoord,:);

resinit=OCMATFTE.bcinitial(depvara,OCMATFTE.initialcoordinate,OCMATFTE.initialstate,modelpar,OCMATCONT.HE.arcarg(1));
resinit=[resinit;depvarafund(:,1)-OCMATFTE.fundsolinit];
restrans=OCMATFTE.bctransversality(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end));
if ~isempty(OCMATFTE.fixendstatecoord)
    restrans(OCMATFTE.fixendstatecoord,1)=depvarb(OCMATFTE.fixendstatecoord,end)-OCMATFTE.endstate;
end
if OCMATFTE.optimalhorizon
    restrans=[restrans; ...
        OCMATFTE.bcoptimalhorizon(timehorizon,depvarb,modelpar,OCMATCONT.HE.arcarg(end))];
end

if OCMATFTE.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATFTE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATFTE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii);
        depvarbfund(:,ii)-depvarafund(:,ii+1)];
end
res=[resinit;resconnec;restrans];

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
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if ~OCMATFTE.optimalhorizon
        arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
    else
        arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
    end
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
    if isfinite(OCMATFTE.maxhorizon)
        violationmat=OCMATFTE.maxhorizon-arctime(end)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
        if violationmat
            counter=counter+1;
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='maxhorizon';
            infoS(counter).cols=1;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=OCMATFTE.maxhorizon-arctime(end);
            infoS(counter).minval=OCMATFTE.maxhorizon-arctime(end);
            b=min([b infoS(counter).minval]);
        end
    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATFTE OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
figure(1)
h=OCMATFTE.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));

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
[t,y,z,freepar]=drearr(tmesh,coeff);

failed=[];
for ii=id
    switch ii
        case 1
            out=OCMATFTE.targetvalue-freepar(OCMATCONT.HE.numparameter);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATFTE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATFTE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
if ~OCMATFTE.optimalhorizon
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).' OCMATFTE.truncationtime];
else
    arctime=[OCMATFTE.initialtime freepar(OCMATFTE.switchtimecoord).'  freepar(OCMATFTE.optimalhorizoncoord)];
end
out.arcinterval=arctime;
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATFTE.initialtime;
out.timehorizon=arctime(end);
out.modelparameter=OCMATFTE.parametervalue;
out.modelparameter(OCMATFTE.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremalpF4ft';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATFTE.switchtimecoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.solverinfo.varyparameterindex=OCMATFTE.varyparameterindex;
out.solverinfo.initialcoordinate=OCMATFTE.initialcoordinate;
out.solverinfo.fixendstatecoord=OCMATFTE.fixendstatecoord;
out.solverinfo.objectivevaluecalc=OCMATFTE.objectivevaluecalc;
out.solverinfo.optimalhorizon=OCMATFTE.optimalhorizon;

% fundamental solution specific variables
out.solverinfo.fundcoord=OCMATFTE.fundcoord;
out.solverinfo.fundmatcoord=OCMATFTE.fundmatcoord;
out.solverinfo.fundnumcoord=OCMATFTE.fundnumcoord;

if OCMATFTE.optimalhorizon
    out.solverinfo.optimalhorizoncoord=OCMATFTE.optimalhorizoncoord;
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

modelpar=OCMATFTE.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        % handle case of pure state constraints
    otherwise
end
modelpar(OCMATFTE.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);


%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATFTE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATFTE=OCMATFTE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATFTE.basicglobalvarfilename '4extremalpF4ft'],'MODELINFO')
    end
    save([OCMATFTE.basicresultfilename '4extremalpF4ft'],'sout','bvpout')
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