function out=extremal2ftae()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;

out{7}=@defaultprocessor;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
out{11}=@operatorpfrechet;
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{17}=@adapt;
out{18}=@dataadaptation;
out{19}=@workspaceadapt;
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

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCMATAE.freehorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freehorizoncoord)];
    %arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.switchtimecoord(end))+OCMATAE.freehorizon];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATAE.objectivevaluecalc
    dxdt(OCMATAE.objectivevaluecoord,:)=dtds*OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
end
%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCMATAE.freehorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freehorizoncoord)];
    %arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.switchtimecoord(end))+OCMATAE.freehorizon];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
if OCMATAE.inftimetransformation && arc==OCMATCONT.HE.numarc
    t=arctime(arc)+log(arc-s)/OCMATAE.inftimetransformation;
    dtds=-1./(OCMATAE.inftimetransformation*(arc-s));
else
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    dtds=diffarctime(arc);
end
J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATAE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATAE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(1).numeq,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATAE.autonomous
        Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if OCMATAE.objectivevaluecalc
        dxdt(OCMATAE.objectivevaluecoord,:)=OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
        Jt(OCMATAE.objectivevaluecoord,:)=OCMATAE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg);
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATCONT.DOMAINDDATA(1).eqcoord,OCMATAE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(1).eqcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        if OCMATAE.inftimetransformation
            Jpar(OCMATCONT.DOMAINDDATA(1).eqcoord,arc-1)=0;
            %Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=-1./((exp(OCMATAE.inftimetransformation*arctime(arc))*(arc-s)))*dxdt;
        else
            Jpar(OCMATCONT.DOMAINDDATA(1).eqcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    end
end
if OCMATAE.freehorizon
    dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
    Jpar(:,OCMATAE.freehorizoncoord)=dxdt;
end
%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP

switchtimes=freepar(OCMATAE.switchtimecoord);
resconnec=[];

initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
resinit=OCMATAE.bcinitial(depvara,OCMATAE.targetcoordinate,initialstate,modelpar,OCMATCONT.HE.arcarg(1));
resasym=OCMATAE.bcasymptotic(depvarb,[],OCMATAE.saddlepoint);
if OCMATAE.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    if OCMATAE.objectivevaluecalc
        resconnec=[resconnec;depvarb(end,ii)-depvara(end,ii+1)];
    end
end
res=[resinit;resconnec;resasym];


%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATAE
[Ja Jb Jpar]=OCMATAE.bcjacobian(depvara,depvarb,freepar,modelpar);

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
figure(1)
h=OCMATAE.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATAE
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
if OCMATAE.stopcriterion
    fprintf(1,' Tangent direction parameter: %g\n',tangent(end));
end
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
global OCMATCONT OCMATAE
failed=[];
out=[];
if isempty(coeff)
    failed=1;
    return
end

for ii=id
    switch ii
        case 1
            if ~OCMATAE.stopcriterion
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
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
global OCMATCONT OCMATAE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if ~OCMATAE.freehorizon
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freehorizoncoord)];
    out.solverinfo.freehorizoncoord=OCMATAE.freehorizoncoord;
    %arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.switchtimecoord(end))+OCMATAE.freehorizon];
end
if ~OCMATAE.freehorizon
    out.timehorizon=OCMATAE.truncationtime;
else
    out.timehorizon=freepar(OCMATAE.freehorizoncoord);
end

out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.conttype='extremal2ftae';

% add solver method specific information to make it consistent with MATLAB
% syntax
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATAE.PD_phi=p'/norm(p);
        OCMATAE.PD_psi=Q(:,end);
        s.data.phi=OCMATAE.PD_phi(:);
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


%-----------------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATAE

modelpar=OCMATAE.parametervalue;
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
%[ttmesh,y,z,freepar,modelpar]=OCMATAE.drearr(tmesh,y,z,freepar,modelpar);
%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATAE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATAE=OCMATAE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATAE.basicglobalvarfilename '4extremal2ftae'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremal2ftae'],'sout','bvpout')
catch
    failed=1;
end

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if ~OCMATAE.freehorizon
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.freehorizoncoord)];
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATAE.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATAE.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATAE

discretizationdata=OCMATAE.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATAE

pathname=OCMATAE.datapath();

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATAE
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
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
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end

