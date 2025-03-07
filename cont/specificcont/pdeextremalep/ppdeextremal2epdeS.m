function out=ppdeextremal2epdeS()
% for a simplified continuation 
out{1}=@operatoreq;
out{2}=@frechetder;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{5}{1}=@odejac;
out{5}{2}=[];%@bcjac;
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
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
out{25}=@deval;
out{26}=@saveintermediate;
out{27}=@datapath;
out{28}=@domaindiscretization;
out{30}=@printcontinuation;

function res=operatoreq(tmesh,coeff,tangent,odefun,bcfun,icfun)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
res=calcpde_RHS(t,y,freepar,modelpar,odefun,bcfun);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac)

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
J=calcpde_RHSJac(t,y,freepar,modelpar,odefun,bcfun,odejacfun,bcjacfun);

% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end
%J=Jnum;

function [opt out]=probleminit(tmesh,coeff,tangent)
global OCMATCONT OCMATPPDESD
opt=setocoptions('OCCONTARG','WorkSpace',1);
opt=setocoptions(opt,'OCCONTARG','Locators',[1 0]);
out=0;


function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,freepar)
global OCMATPPDESD OCMATCONT OCBVP

modelpar=OCMATPPDESD.parametervalue;
x=OCMATPPDESD.points;
if ~OCMATPPDESD.movinghorizon
    arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' OCMATPPDESD.truncationtime];
else
    arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' freepar(OCMATPPDESD.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcid=0;
arc=1;
t=diffarctime(arc)*s+(arctime(arc+1)-diffarctime(arc)*arc);
dtds=diffarctime(arc);
dxdt=zeros(size(depvar));
for ii=1:length(s)
    depvarii=depvar(OCMATPPDESD.totalcoordinate,ii);
    F=OCMATPPDESD.canonicalsystem(t,x,depvarii(OCMATPPDESD.coeffidx),modelpar,arcid);
    %dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(OCMATPPDESD.femop.M*F-(OCMATPPDESD.femop.K-OCMATPPDESD.femop.Kadv)*depvar(OCMATPPDESD.totalcoordinate,ii));
    %dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(OCMATPPDESD.femop.invM*F-OCMATPPDESD.femop.invMKmKadv*depvar(OCMATPPDESD.totalcoordinate,ii));
    if ~OCMATPPDESD.MassMatrix
        dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(F-OCMATPPDESD.femop.invMKmKadv*depvar(OCMATPPDESD.totalcoordinate,ii));
    else
        dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(OCMATPPDESD.femop.M*F-OCMATPPDESD.femop.K*depvar(OCMATPPDESD.totalcoordinate,ii));
    end
    if OCMATPPDESD.objectivevaluecalc
        depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,depvar(OCMATPPDESD.totalcoordinate,ii));
        dxdt(OCMATPPDESD.objectivevaluecoord,ii)=dtds*sum(OCMATPPDESD.trianglearea.*OCMATPPDESD.objectivefunction(t(ii),OCMATPPDESD.points,depvarint,modelpar,arcid));
    end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jmodelpar]=odejac(s,depvar,freepar)
global OCMATPPDESD OCMATCONT OCBVP

modelpar=OCMATPPDESD.parametervalue;
x=OCMATPPDESD.points;
if ~OCMATPPDESD.movinghorizon
    arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' OCMATPPDESD.truncationtime];
else
    arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' freepar(OCMATPPDESD.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcid=0;
arc=1;
t=diffarctime(arc)*s+(arctime(arc+1)-diffarctime(arc)*arc);
dtds=diffarctime(arc);

Fu=OCMATPPDESD.canonicalsystemjacobian(t,x,depvar(OCMATPPDESD.coeffidx),modelpar,arcid);
if ~OCMATPPDESD.MassMatrix
    J=dtds*(Fu-OCMATPPDESD.femop.inexactinvMKmKadv);
else
    J=dtds*(OCMATPPDESD.femop.M*Fu-OCMATPPDESD.femop.K);
end

if OCMATPPDESD.objectivevaluecalc
    depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,depvar(OCMATPPDESD.totalcoordinate));
    OJ=(dtds.*OCMATPPDESD.triangleareajac.*OCMATPPDESD.objectivefunctionjacobian(t,x,depvarint,modelpar,arcid)*OCMATPPDESD.JacInt.').';
    OJ=[OJ(:).' 0];
    J=[J zeros(OCBVP.numode-1,1); ...
        OJ];
end


Jmodelpar=zeros(OCBVP.numode,OCBVP.numparameter);
if OCMATPPDESD.movinghorizon
    F=OCMATPPDESD.canonicalsystem(t,x,depvar(OCMATPPDESD.coeffidx),modelpar,arcid);
    if ~OCMATPPDESD.MassMatrix
        dxdt(OCMATPPDESD.totalcoordinate)=(F-OCMATPPDESD.femop.invMKmKadv*depvar(OCMATPPDESD.totalcoordinate));
    else
        dxdt(OCMATPPDESD.totalcoordinate)=(OCMATPPDESD.femop.M*F-OCMATPPDESD.femop.K*depvar(OCMATPPDESD.totalcoordinate));
    end
    if OCMATPPDESD.objectivevaluecalc
        depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,depvar(OCMATPPDESD.totalcoordinate));
        dxdt(OCMATPPDESD.objectivevaluecoord,ii)=sum(OCMATPPDESD.trianglearea.*OCMATPPDESD.objectivefunction(t,OCMATPPDESD.points,depvarint,modelpar,arcid));
    end
    Jmodelpar(:,1)=dxdt;
end
% options.diffvar=2; 
% options.vectvars=[];
% dFdy=numjaccsd(@ode,{s,depvar,freepar},size(J,1),options);
% options.diffvar=3; 
% dFdpar=numjaccsd(@ode,{s,depvar,freepar},size(Jmodelpar,1),options);
% J-dFdy
% 
%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar)
global OCMATPPDESD OCMATCONT OCBVP

initialstate=OCMATPPDESD.startdistribution+freepar(end)*OCMATPPDESD.continuationvector;

rescontpar=sum(([depvara;freepar]-OCMATCONT.LastIntialDistribution-OCMATCONT.InitialSecant).*OCMATCONT.InitialSecant);
resinit=OCMATPPDESD.bcinitial(depvara(OCMATPPDESD.totalcoordinate),OCMATPPDESD.targetcoordinate,initialstate);
resasym=OCMATPPDESD.bcasymptotic(depvarb(OCMATPPDESD.totalcoordinate),OCMATPPDESD.asymptoticmatrix,OCMATPPDESD.saddlepoint);
if OCMATPPDESD.movinghorizon
    resasym=[resasym; ...
        sqrt(sum((OCMATPPDESD.saddlepoint-depvarb(OCMATPPDESD.totalcoordinate,end)).^2))-OCMATPPDESD.distance];
end
if OCMATPPDESD.objectivevaluecalc
    resinit=[resinit;depvara(OCMATPPDESD.objectivevaluecoord,1)];
end

res=[resinit;resasym;rescontpar];

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATPPDESD OCBVP

b=0;
infoS=[];
labelS='';
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff)
global OCMATPPDESD OCMATCONT
[t,y,freepar,modelpar]=drearr(tmesh,coeff); 
%sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATPPDESD.plotcontinuation(t,y,modelpar,freepar);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)

idx=[];
if isempty(coeff)
    return
end
[t,y,freepar]=drearr(tmesh,coeff); 
fprintf(1,' Continuation parameter: %g\n',freepar(end));

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATPPDESD
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
global OCMATCONT OCMATPPDESD
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
global OCMATCONT OCMATPPDESD OCBVP

[t,y,freepar,modelpar]=drearr(tmesh,coeff); 

out=transform2nativematlab(tmesh,coeff);
out.t=t;
out.freeparameter=freepar;
if ~OCMATPPDESD.movinghorizon
    timeinterval=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' OCMATPPDESD.truncationtime];
else
    timeinterval=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' freepar(OCMATPPDESD.movinghorizoncoord)];
end
out.discretizationinfo.timeinterval=timeinterval;
out.optidentifier=[];
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;
out.continuationparameter=freepar(end);

out.discretizationinfo.method=OCMATCONT.bvpmethod;
out.discretizationinfo.coeff=coeff;
out.discretizationinfo.tangent=tangent;
out.discretizationinfo.tmesh=tmesh;
out.discretizationinfo.parameters=freepar;
out.discretizationinfo.coord.contparameter=OCMATCONT.HE.contparametercoord;
out.discretizationinfo.coord.switchtime=OCMATPPDESD.switchtimecoord;
out.discretizationinfo.coord.totaldepvar=OCMATPPDESD.totalcoordinate;
out.discretizationinfo.coord.meshidx=OCMATPPDESD.coeffidx;
if OCMATPPDESD.objectivevaluecalc
    out.discretizationinfo.coord.objectivevalue=OCMATPPDESD.objectivevaluecoord;
else
    out.discretizationinfo.coord.objectivevalue=[];
end
if OCMATPPDESD.movinghorizon
    out.discretizationinfo.coord.movinghorizoncoord=OCMATPPDESD.movinghorizoncoord;
    out.discretizationinfo.distance=OCMATPPDESD.distance;
end
out.discretizationinfo.pathtype=OCMATPPDESD.pathtype;
out.discretizationinfo.points=OCMATPPDESD.points;
out.discretizationinfo.edges=OCMATPPDESD.edges;
out.discretizationinfo.triangles=OCMATPPDESD.triangles;
out.spacegeometry.csgregion=OCMATPPDESD.geo;
out.spacegeometry.csgbooleantable=OCMATPPDESD.bt;
out.boundarycondition=OCMATPPDESD.boundarycondition;

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATPPDESD
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATPPDESD.PD_phi=p'/norm(p);
        OCMATPPDESD.PD_psi=Q(:,end);
        s.data.phi=OCMATPPDESD.PD_phi(:);
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
function [tmesh,y,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATPPDESD OCBVP

modelpar=OCMATPPDESD.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        N=length(tmesh);
        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:N*(OCBVP.numode),OCBVP.numode,N);
        OCMATCONT.HE.parametercoord=N*OCBVP.numode+(1:OCBVP.numparameter).';
        OCMATCONT.HE.contparametercoord=OCMATCONT.HE.parametercoord(end);
        OCMATCONT.LastIntialDistributionIndex=[OCMATCONT.HE.DDATA.meshvalcoord(:,1);OCMATCONT.HE.parametercoord(:)];
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    case {'mtom0','tom'}
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
global OCMATCONT OCMATPPDESD OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATPPDESD=OCMATPPDESD;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATPPDESD.basicglobalvarfilename '4ppdeextremal2epdeS'],'MODELINFO')
    end
    save([OCMATPPDESD.basicresultfilename '4ppdeextremal2epdeS'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATPPDESD

pathname=OCMATPPDESD.datapath();

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATPPDESD

discretizationdata=OCMATPPDESD.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

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
        paramcoeff=n+(1:OCBVP.numparameter);
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
    case {'mtom0','tom'}
        sol.y=y;
        sol.x=t;
        sol.parameters=freepar;
        sol.solver=OCMATCONT.bvpmethod;
    otherwise
        return
end

%-----------------------------------------------------------------
function [dFdy,nfevals,nfcalls]=numjaccsd(fun,Fargs,nF,options)
% NUMJACCSD    Complex Step Jacobian
% based on 'jacobiancsd' by Yi Cao at Cranfield University, 02/01/2008 and
% uses the structure of odenumjac.
%
% J = NUMJACCSD(F,FARGS,NF,OPTIONS) returns the numerical (NF x N) Jacobian
% matrix of a NF-vector function, F(FARGS{:}) at the reference point, 
% Y=FARGS{OPTIONS.DIFFVAR} (N-vector). 
%
% The structure OPTIONS must have the following fields: DIFFVAR, VECTVARS.
% The field OPTIONS.DIFFVAR is the index of the  differentiation variable,
% Y = FARGS{DIFFVAR}. For a function F(t,x), set DIFFVAR to 1 to compute
% DF/Dt, or to 2 to compute DF/Dx. Set OPTIONS.VECTVAR to the indices of
% vectorized arguments: VECTVAR = [2] indicates that  F(t,[x1 y2 ...])
% returns [F(t,x1) F(t,x2) ...], while VECTVAR = [1,2] indicates that F([t1
% t2 ...],[x1 x2 ...]) returns [F(t1,x1) F(t2,x2) ...].
%   
% [DFDY,NFEVALS,NFCALLS] = ODENUMJAC(...) returns the number of values
% F(FARGS{:}) computed while forming dFdY (NFEVALS) and the number of calls
% to the function F (NFCALLS). If F is not vectorized, NFCALLS equals
% NFEVALS.  
% Dieter Grass

% Options
diffvar = options.diffvar; 
vectvar = options.vectvars;

% The differentiation variable.
y  = Fargs{diffvar};
classY = class(y);

ny = length(y);

dFdy=zeros(nF,ny);                   % allocate memory for the Jacobian matrix
del = (y + ny*eps(classY)*i) - y;
h=imag(del);
ydel = y(:,ones(1,ny)) + diag(del);
if isempty(vectvar)
    for ii=1:ny                      % loop for each independent variable
        dFdy(:,ii)=imag(fun(Fargs{1:diffvar-1},ydel(:,ii),Fargs{diffvar+1:end}))/h(ii);     % complex step differentiation
    end
    nfcalls = ny;                       % stats
else
    Fargs_expanded = Fargs;
    Fargs_expanded{diffvar} = ydel;
    vectvar = setdiff(vectvar,diffvar);
    for ii=1:length(vectvar)
      Fargs_expanded{vectvar(ii)} = repmat(Fargs{vectvar(ii)},1,ny);
    end
    dFdy=imag(fun(Fargs_expanded{:}))./h(:,ones(1,ny)); 
    nfcalls = 1;                       % stats
end
nfevals = ny;                         % stats (at least one per loop)
