function out=ppdeextremalt2epdeS()

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
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATPPDESD OCMATCONT OCBVP

x=OCMATPPDESD.points;
arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' freepar(OCMATCONT.HE.numparameter)];
diffarctime=diff(arctime);
arcid=0;

t=diffarctime(arc)*s+(arctime(arc+1)-diffarctime(arc)*arc);
dtds=diffarctime(arc);
dxdt=zeros(size(depvar));
for ii=1:length(s)
    depvarii=depvar(OCMATPPDESD.totalcoordinate,ii);
    F=OCMATPPDESD.canonicalsystem(t,x,depvarii(OCMATPPDESD.coeffidx),modelpar,arcid);
    %dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(OCMATPPDESD.femop.M*F-(OCMATPPDESD.femop.K-OCMATPPDESD.femop.Kadv)*depvar(OCMATPPDESD.totalcoordinate,ii));
    %dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(OCMATPPDESD.femop.invM*F-OCMATPPDESD.femop.invMKmKadv*depvar(OCMATPPDESD.totalcoordinate,ii));
    dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(F-OCMATPPDESD.femop.invMKmKadv*depvar(OCMATPPDESD.totalcoordinate,ii));
    %dxdt(OCMATPPDESD.totalcoordinate,ii)=dtds*(OCMATPPDESD.femop.M\F-OCMATPPDESD.femop.invMKmKadv*depvar(OCMATPPDESD.totalcoordinate,ii));
    if OCMATPPDESD.objectivevaluecalc
        depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,depvar(OCMATPPDESD.totalcoordinate,ii));
        dxdt(OCMATPPDESD.objectivevaluecoord,ii)=dtds*sum(OCMATPPDESD.trianglearea.*OCMATPPDESD.objectivefunction(t(ii),OCMATPPDESD.points,depvarint,modelpar,arcid));
    end
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jmodelpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATPPDESD OCMATCONT OCBVP
x=OCMATPPDESD.points;

arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' freepar(OCMATCONT.HE.numparameter)];
diffarctime=diff(arctime);
arcid=0;
t=diffarctime(arc)*s+(arctime(arc+1)-diffarctime(arc)*arc);
dtds=diffarctime(arc);

Fu=OCMATPPDESD.canonicalsystemjacobian(t,x,depvar(OCMATPPDESD.coeffidx),modelpar,arcid);
%J=dtds*(OCMATPPDESD.femop.M*Fu-(OCMATPPDESD.femop.K-OCMATPPDESD.femop.Kadv));
%J=dtds*(OCMATPPDESD.femop.invM*Fu-OCMATPPDESD.femop.invMKmKadv);
J=dtds*(Fu-OCMATPPDESD.femop.invMKmKadv);
%J=dtds*(OCMATPPDESD.femop.M\Fu-OCMATPPDESD.femop.invMKmKadv);
if OCMATPPDESD.objectivevaluecalc
    %zerov=zeros(OCMATPPDESD.statecostatenumcoordinate(end),1);
    depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,depvar(OCMATPPDESD.totalcoordinate));
    OJ=(dtds.*OCMATPPDESD.triangleareajac.*OCMATPPDESD.objectivefunctionjacobian(t,x,depvarint,modelpar,arcid)*OCMATPPDESD.JacInt.').';
    OJ=[OJ(:).' 0];
    J=[J zeros(OCBVP.neqn-1,1); ...
        OJ];
    numJacOpt.diffvar=2;
%     numJacOpt.vectvars=[];
%     Jnum=numjaccsd(@calcobjective,{s,depvar,arc,freepar,modelpar},1,numJacOpt);
%     if any(abs(Jnum-OJ(:).')>1e-1)
%         Jnum-OJ(:).'
%     end
%     J=[J; ...
%         OJ(:).'];
%     J=[J zeros(OCBVP.neqn-1,1); ...
%         Jnum];
end


Jmodelpar=zeros(OCBVP.neqn,OCMATCONT.HE.numparameter);



function out=pdeboundarycondition()
global OCMATPPDESD
out=OCMATPPDESD.boundarycondition;

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATPPDESD OCMATCONT OCBVP
resconnec=[];

initialstate=OCMATPPDESD.startdistribution;
resinit=OCMATPPDESD.bcinitial(depvara(OCMATPPDESD.initialcoordinate),OCMATPPDESD.initialcoordinate,initialstate);
resasym=OCMATPPDESD.bcasymptotic(depvarb(OCMATPPDESD.totalcoordinate),OCMATPPDESD.asymptoticmatrix,OCMATPPDESD.saddlepoint);
if OCMATPPDESD.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end

res=[resinit;resconnec;resasym];

%----------------------------------------------------------------
function [res,varargout]=residual(tmesh,coeff,rhs,odefun)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        res=residualpde_bvp6c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
    case 'bvp4c'
        res=residualpde_bvp4c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
end

%----------------------------------------------------------------
function res=maxresidual(tmesh,coeff)
global OCMATCONT

switch OCMATCONT.bvpmethod
    case 'bvp5c'
        [t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
        res=maxresidual_bvp5c(t,y);
    otherwise
        res=0;
end

%----------------------------------------------------------------
function b=verifyresidual(maxres)
global OCMATCONT

switch OCMATCONT.bvpmethod
    case 'bvp5c'
        b=verifyresidual_bvp5c(maxres);
    otherwise
        b=maxres < OCMATCONT.OPTIONS.meshadaptreltol;
end

%----------------------------------------------------------------
function flag=prepare4meshadapt(tmesh,coeff)
global OCMATCONT

switch OCMATCONT.bvpmethod
    case 'bvp5c'
        [t,y]=drearr(tmesh,coeff);
        prepare4meshadapt_bvp5c(t,y);
        flag=1;
    otherwise
        flag=1;
end

%----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=meshadaptation(tmesh,coeff,tangent,res,canRemovePoints,ode)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp6c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp5c(t,y,tangent,res,canRemovePoints,freepar,modelpar,ode);
        coeffnew=[ynew(:)];
    case 'bvp4c'
        [tmeshnew,ynew,tangentnew]=meshadaptation_bvp4c(t,y,tangent,res,canRemovePoints);
        coeffnew=[ynew(:);freepar(:)];
end

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATPPDESD OCBVP

b=0;
infoS=[];
labelS='';
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATPPDESD OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
%sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATPPDESD.plotcontinuation(t,y,modelpar,freepar);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATPPDESD
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATPPDESD.targettype
    case 'd'
        fprintf(1,' Distance |yhat-yend|: %3.7g\n',norm(OCMATPPDESD.saddlepoint-y(:,end))-OCMATPPDESD.targetvalue);
    case 'T'
        fprintf(1,' Difference Time horizon: %g\n',freepar(OCMATPPDESD.HE.numparameter)-OCMATPPDESD.targetvalue);
    otherwise
        out=[];
end


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
global OCMATCONT OCMATPPDESD

failed=[];
for ii=id
    switch ii
        case 1
            switch OCMATPPDESD.targettype
                case 'd'
                    out=norm(OCMATPPDESD.saddlepoint-y(:,end))-OCMATPPDESD.targetvalue;
                case 'T'
                    out=freepar(OCMATCONT.HE.numparameter)-OCMATPPDESD.targetvalue;
                otherwise
                    out=[];
            end
        otherwise
            out=[];
            failed=1;
    end
end
%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATPPDESD OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 

out=transform2nativematlab(tmesh,coeff,OCMATPPDESD.parametervalue);
out.t=out.x;
out.freeparameter=freepar;
out.discretizationinfo.timeinterval=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' freepar(OCMATCONT.HE.numparameter)];
out.optidentifier=[];
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.discretizationinfo.method=OCMATCONT.bvpmethod;
out.discretizationinfo.coeff=coeff;
out.discretizationinfo.tangent=tangent;
out.discretizationinfo.tmesh=tmesh;
out.discretizationinfo.parameters=freepar;
out.discretizationinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.discretizationinfo.coord.contparameter=OCMATCONT.HE.contparametercoord;
out.discretizationinfo.coord.switchtime=OCMATPPDESD.switchtimecoord;
out.discretizationinfo.coord.totaldepvar=OCMATPPDESD.totalcoordinate;
out.discretizationinfo.coord.meshidx=OCMATPPDESD.coeffidx;
if OCMATPPDESD.objectivevaluecalc
    out.discretizationinfo.coord.objectivevalue=OCMATPPDESD.objectivevaluecoord;
else
    out.discretizationinfo.coord.objectivevalue=[];
end
out.discretizationinfo.pathtype=OCMATPPDESD.pathtype;
out.discretizationinfo.points=OCMATPPDESD.points;
out.discretizationinfo.edges=OCMATPPDESD.edges;
out.discretizationinfo.triangles=OCMATPPDESD.triangles;
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.discretizationinfo.yp=out.yp;
        out.discretizationinfo.ypmid=out.ypmid;
        out=rmfield(out,{'yp','ypmid'});
    case 'bvp5c'
    case 'bvp4c'
        out.discretizationinfo.yp=out.yp;
        out=rmfield(out,'yp');
end
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
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATPPDESD
z=[];
modelpar=OCMATPPDESD.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
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
        save([OCMATPPDESD.basicglobalvarfilename '4ppdeextremal2epde'],'MODELINFO')
    end
    save([OCMATPPDESD.basicresultfilename '4ppdeextremal2epde'],'sout','bvpout')
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

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATPPDESD
switch OCMATCONT.bvpmethod
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATPPDESD.parametervalue);
        nN=n*N;
        % under the giving assumptions the solution is continuous even at
        % the border of the arcs. Therefore, the correspondent warning of
        % deval is turned 'off'
        warning('off','MATLAB:deval:NonuniqueSolution');
        %sol.yp=sol.discretizationinfo.yp;
        %try
        %    sol.ypmid=sol.discretizationinfo.ypmid;
        %end
        ynew=devalbvpoc(sol,tmeshnew);
        warning('on','MATLAB:deval:NonuniqueSolution');
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,n,N,nN);

        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATPPDESD.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.discretizationinfo.yp;
        sol.idata.ymid=sol.discretizationinfo.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end

function [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,modelpar)
global OCMATCONT OCBVP
N=numel(tmesh);
sol=[];
freepar=[];
n=[];
switch OCMATCONT.bvpmethod
    case {'bvp4c','bvp6c'}
        n=OCBVP.numode;
        nN=n*N;
        y=reshape(coeff(1:nN),n,N);
        freepar=coeff(nN+(1:OCMATCONT.HE.numparameter));
        if strcmp(OCMATCONT.bvpmethod,'bvp4c')
            yp=calcdata_bvp4c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
        elseif strcmp(OCMATCONT.bvpmethod,'bvp6c')
            [yp,ypmid]=calcdata_bvp6c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
            sol.ypmid=ypmid;
        end
        sol.x=tmesh;
        sol.y=y;
        sol.yp=yp;
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
    otherwise
        return
end
if OCMATCONT.HE.numparametermc>0
    sol.parameters=freepar(1:OCMATCONT.HE.numparametermc);
end

%-------------------------------------------------------------------------
function o=calcobjective(s,depvar,arc,freepar,modelpar)
global OCMATPPDESD

arctime=[OCMATPPDESD.initialtime freepar(OCMATPPDESD.switchtimecoord).' OCMATPPDESD.truncationtime];
diffarctime=diff(arctime);
arcid=0;

t=diffarctime(arc)*s+(arctime(arc+1)-diffarctime(arc)*arc);
dtds=diffarctime(arc);
o=zeros(1,length(s));
for ii=1:length(s)
    if OCMATPPDESD.objectivevaluecalc
        depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,depvar(OCMATPPDESD.totalcoordinate,ii));
        o(ii)=dtds*sum(OCMATPPDESD.trianglearea.*OCMATPPDESD.objectivefunction(t(ii),OCMATPPDESD.points,depvarint,modelpar,arcid));
    end
end
