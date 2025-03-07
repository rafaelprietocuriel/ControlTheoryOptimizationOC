function out=extremal2epS()

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
out{12}{1}=@residual;
out{12}{2}=@maxresidual;
out{12}{3}=@verifyresidual;
out{12}{4}=@prepare4meshadapt;
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{17}=@adapt;
out{18}=@meshadaptation;
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
function dxdt=ode(s,depvar,freepar)
global OCMATAE OCMATCONT OCBVP
modelpar=OCMATAE.parametervalue;
arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
arc=1;
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
function [J Jpar]=odejac(s,depvar,freepar)
global OCMATAE OCMATCONT OCBVP
modelpar=OCMATAE.parametervalue;
arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
diffarctime=diff(arctime);
arc=1;
arcindex=arc;
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);

J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATAE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATAE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCBVP.numparameter);
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
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    else
        if OCMATAE.inftimetransformation
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=0;
            %Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=-1./((exp(OCMATAE.inftimetransformation*arctime(arc))*(arc-s)))*dxdt;
        else
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
        end
    end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [H Hpar]=odehess(s,depvar,arc,freepar)
global OCMATAE OCMATCONT OCBVP
modelpar=OCMATAE.parametervalue;
if isempty(OCMATAE.movinghorizon)
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.switchtimecoord(end))+OCMATAE.movinghorizon];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
H=dtds*OCMATAE.canonicalsystemhessian(t,depvar,modelpar,arcarg);
%H(end+(1:OCBVP.numparameter),:,:)=0;
Hpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.DOMAINDDATA(arcindex).numeq,OCBVP.numparameter);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar)
global OCMATAE OCMATCONT OCBVP
modelpar=OCMATAE.parametervalue;
switchtimes=freepar(OCMATAE.switchtimecoord);
if OCMATAE.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATAE.jumpcostateindex)=freepar(OCMATAE.jumpcostatecoord);
end
resconnec=[];

initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
%rescontpar=sum(([depvara(OCMATAE.targetcoordinate,1);freepar]-OCMATCONT.LastIntialDistribution([OCMATAE.targetcoordinate end])-OCMATCONT.InitialSecant([OCMATAE.targetcoordinate end])).*OCMATCONT.InitialSecant([OCMATAE.targetcoordinate end]));
rescontpar=sum(([depvara;freepar]-OCMATCONT.LastIntialDistribution-OCMATCONT.InitialSecant).*OCMATCONT.InitialSecant);
resinit=OCMATAE.bcinitial(depvara,OCMATAE.targetcoordinate,initialstate,modelpar,OCMATCONT.HE.arcarg(1));
resasym=OCMATAE.bcasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint);
if OCMATAE.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end

if OCMATAE.stateconstraint
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,jump,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
else
    for ii=1:numel(OCMATCONT.HE.arcarg)-1
        resconnec=[resconnec;
            OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
            OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
    end
end
res=[resinit;resconnec;resasym;rescontpar];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar)
global OCMATAE OCMATCONT OCBVP
switchtimes=freepar(OCMATAE.switchtimecoord);
modelpar=OCMATAE.parametervalue;
if OCMATAE.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATAE.jumpcostateindex)=freepar(OCMATAE.jumpcostatecoord);
end
domainddata=OCMATCONT.DOMAINDDATA;

Jpar=zeros(OCMATCONT.HE.totalnumboundarycondition,OCBVP.numparameter);
Ja=zeros(OCMATCONT.HE.totalnumboundarycondition-OCBVP.numparameter+OCMATCONT.codimension);
Jb=Ja;
colcounter=0;
rowcounter=0;
arcindex=OCMATCONT.HE.arcindex(1);
colidx_start=colcounter+1;
colcounter=colcounter+domainddata(arcindex).numeq;
[Japart Jbpart Jpartmp]=OCMATAE.bcjacobianinitial(depvara,freepar,OCMATCONT.HE.arcarg(1),OCMATAE.targetcoordinate,OCMATAE.continuationvector);
if OCMATAE.objectivevaluecalc
    Japart(OCMATAE.objectivevaluecoord,OCMATAE.objectivevaluecoord)=1;
end

Jpar(1:OCMATCONT.HE.numinitialcondition,:)=Jpartmp;
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numinitialcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;
for arc=1:OCMATCONT.HE.numarc-1
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    Jaresetpart=OCMATAE.jacobianreset(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,-1);
    Jbresetpart=OCMATAE.jacobianreset(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,0);
    Jaguardpart=OCMATAE.jacobianguard(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,-1);
    Jbguardpart=OCMATAE.jacobianguard(depvara,depvarb,modelpar,switchtimes,arcarg,OCMATCONT.HE.edge,arc,0);
    
    % Ja
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+domainddata(arcindex).numode;
    Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbresetpart;
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+1;
    Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbguardpart;
    
    colidx_start=colcounter+1;
    colcounter=colcounter+domainddata(arcindex).numeq;
    % Jb
    rowcounter=rowcounter-domainddata(arcindex).numode-1; % 
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+domainddata(arcindex).numode;
    Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Jaresetpart;
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+1;
    Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Jaguardpart;
end

[Japart Jbpart  Jpartmp]=OCMATAE.bcjacobianasymptotic(depvarb,OCMATCONT.HE.arcarg(OCMATCONT.HE.numarc),OCMATAE.asymptoticmatrix',OCMATAE.saddlepoint);
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numendcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;


%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if isempty(OCMATAE.movinghorizon)
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.switchtimecoord(end))+OCMATAE.movinghorizon];
    end
    %arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
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
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff)
global OCMATAE OCMATCONT
[t,y,freepar,modelpar]=drearr(tmesh,coeff);
%sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
figure(1)
h=OCMATAE.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,[]);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,freepar]=drearr(tmesh,coeff); 
fprintf(1,' Continuation parameter: %g\n',freepar(end));

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
            out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATAE OCBVP

[t,y,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff);
out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
%out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
if OCMATAE.inftimetransformation
    out.timehorizon=inf;
else
    if isempty(OCMATAE.movinghorizon)
        out.timehorizon=OCMATAE.truncationtime;
    else
        out.timehorizon=freepar(OCMATAE.switchtimecoord(end))+OCMATAE.movinghorizon;
    end
end
out.modelparameter=OCMATAE.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremal2epS';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATAE.inftimetransformation;
out.solverinfo.pathtype=OCMATAE.pathtype;
%out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATAE.switchtimecoord;
if OCMATAE.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATAE.objectivevaluecoord;
end

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        % reduce Jacobian to size without the continuation parameter
        J(:,OCMATCONT.HE.contparametercoord)=[];
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATAE.PD_phi=p'/norm(p); % right eigenvector
        OCMATAE.PD_psi=Q(:,end); % left eigenvector
        s.data.phi=OCMATAE.PD_phi(:);
        s.data.psi=OCMATAE.PD_psi(:);
        s.data.DFDX=J;
        s.data.laecoefficient=nf_LAE(tmesh,coeff,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac,@odehess,[],s.data.psi,s.data.phi); % quadratic coefficient of center manifold
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
function [tmesh,y,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATAE OCBVP

modelpar=OCMATAE.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        N=length(tmesh);
        OCBVP.N=N;
        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:N*(OCBVP.numode),OCBVP.numode,N);
        OCMATCONT.HE.parametercoord=N*OCBVP.numode+(1:OCBVP.numparameter).';
        OCMATCONT.HE.contparametercoord=OCMATCONT.HE.parametercoord(end);
        OCMATCONT.LastIntialDistributionIndex=[OCMATCONT.HE.DDATA.meshvalcoord(:,1);OCMATCONT.HE.parametercoord(:)];
        ff=find(diff(tmesh)==0);
        OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 ff];
        OCMATCONT.HE.TIMEDDATA.rightarcindex=[ff+1 OCBVP.N];
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        % handle case of pure state constraints
    case 'tom'
        N=length(tmesh);
        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:N*(OCBVP.numode+OCBVP.numparameter),OCBVP.numode+OCBVP.numparameter,N);
        OCMATCONT.HE.DDATA.meshvalcoord(end-OCBVP.numparameter+1:end,:)=[];
        OCMATCONT.HE.parametercoord=OCBVP.numode+(1:OCBVP.numparameter).';
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
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
        save([OCMATAE.basicglobalvarfilename '4extremal2epS'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremal2epS'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATAE

pathname=OCMATAE.datapath();

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATAE

discretizationdata=OCMATAE.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

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
        neqn=n+OCBVP.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end

function c=nf_LAE(tmesh,coeff,ode,bc,odejac,bcjac,odehess,bchess,w,v)
global OCMATCONT OCMATAE
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp4c'
        c=calchnf_bvp4c(t,y,freepar,modelpar,ode,bc,odejac,bcjac,odehess,bchess,w,v);
        %finite difference approximation
%         incr=1e-5;
%         v_y=v(OCMATCONT.HE.DDATA.meshvalcoord);
%         w_y=w(OCMATCONT.HE.DDATA.meshvalcoord);
%         c2=w'*(calc_RHS(t,y-incr*v_y,z,freepar,modelpar,ode,bc,[])+calc_RHS(t,y+incr*v_y,z,freepar,modelpar,ode,bc,[]))/incr^2;
        %         Bvpw=(calc_RHS(t,y-incr*(v_y+w_y),z,freepar,modelpar,ode,bc,[])+calc_RHS(t,y+incr*(v_y+w_y),z,freepar,modelpar,ode,bc,[]))/incr^2;
        %         Bvmw=(calc_RHS(t,y-incr*(v_y-w_y),z,freepar,modelpar,ode,bc,[])+calc_RHS(t,y+incr*(v_y-w_y),z,freepar,modelpar,ode,bc,[]))/incr^2;
        %         Bvw=(Bvpw-Bvmw)/4;
        %[dGy dGpar]=calchessianphi_bvp4c(t,y,freepar,modelpar,ode,bc,odejac,bcjac,odehess,bchess,w,v);
    otherwise
        c=[];
end


function [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff)
global OCMATCONT OCBVP
[t,y,freepar,modelpar]=drearr(tmesh,coeff); 

N=numel(tmesh);
sol=[];
n=[];
switch OCMATCONT.bvpmethod
    case {'bvp4c','bvp6c'}
        n=OCBVP.numode;
        nN=n*N;
        y=reshape(coeff(1:nN),n,N);
        freepar=coeff(nN+(1:OCBVP.numparameter));
        sol.x=tmesh;
        sol.y=y;
        sol.parameters=freepar;
        sol.solver=OCMATCONT.bvpmethod;
    case 'bvp5c'
        n=OCBVP.numode;
        neqn=OCBVP.neqn;
        y=reshape(coeff,neqn,N);
        paramcoeff=n+(1:OCMATCONT.HE.numparameter);
        freepar=coeff(paramcoeff);
        [yp ymid]=calcdata_bvp5c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
        y(paramcoeff,:)=[];
        yp(paramcoeff,:)=[];
        ymid(paramcoeff,:)=[];
        sol.x=tmesh(1:OCBVP.nstages:end);
        sol.y=y(:,1:OCBVP.nstages:end);
        sol.discretizationinfo.yp=yp;
        sol.discretizationinfo.ymid=ymid;
        sol.solver=OCMATCONT.bvpmethod;
    case 'tom'
        sol.y=y;
        sol.x=t;
        sol.parameters=freepar;
        sol.solver=OCMATCONT.bvpmethod;
    otherwise
        return
end
