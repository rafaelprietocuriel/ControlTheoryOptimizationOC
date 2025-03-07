function out=extremaldae4ft()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@dae;
out{4}{2}=@bc;
out{4}{3}=@ic;
out{5}{1}=@daejac;
out{5}{2}=@bcjac;
out{5}{3}=@icjac;
out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@settargetvalue;
out{11}=@operatorpfrechet;
out{12}=@residual;
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{20}=@workspaceadapt;
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
out{25}=@residualdrearr;
out{26}=@saveintermediate;
out{27}=@writelogfile;
out{28}=@datapath;
out{30}=@printcontinuation;

function res=operatoreq(tcolmesh,coeff,tangent,daefun,bcfun,icfun)
global OCMATCONT

[tcolmesh,x,y,z1,freepar,modelpar]=drearr(tcolmesh,coeff,tangent);
% tic
% res0=calcres_dae(tcolmesh,x,y,z1,freepar,modelpar,daefun,bcfun,icfun,tangent);
% t=toc;fprintf('elapsed time: %f\n',t)
tmesh=tcolmesh(1:OCMATCONT.CollocationNumber+1:end);
%tic
res=calcres_dae_new(tmesh,x,y,z1,freepar,modelpar,daefun,bcfun,icfun,tangent);
%t=toc;fprintf('elapsed time: %f\n',t)

% for testing
%testoperatoreq(tcolmesh,coeff,tangent,daefun,bcfun,icfun,res)

function J=frechetder(tcolmesh,coeff,tangent,daefun,bcfun,icfun,daejacfun,bcjacfun,icfunjac)

%tic
[tcolmesh,x,y,z1,freepar,modelpar]=drearr(tcolmesh,coeff,tangent);
J=calcresjac_dae(tcolmesh,x,y,z1,freepar,modelpar,daefun,bcfun,icfun,daejacfun,bcjacfun,icfunjac,tangent);
%t=toc;fprintf('elapsed time: %f\n',t)

% for testing
%testfrechetderivative(tcolmesh,coeff,tangent,daefun,bcfun,icfun,daejacfun,bcjacfun,icfunjac,J,modelpar);

function varargout=probleminit(varargin)
tcolmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tcolmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;

function [F,J]=operatorpfrechet(tcolmesh,coeff,tangent,varargin)
F=operatoreq(tcolmesh,coeff,tangent,varargin{:});
J=frechetder(tcolmesh,coeff,tangent,varargin{:});

function res=residual(tcolmesh,coeff,tangent,tres,daefun)

%tic
[tcolmesh,x,dxdt,freepar,modelpar]=residualrearr(tcolmesh,coeff,tangent,tres);
res=residual_dae(tres,x,dxdt,freepar,modelpar,daefun);

%-------------------------------------------------------------------------
function rhs=dae(tcolmesh,depvar,freepar,modelpar,varargin)
global OCMATFTE
if strcmp(OCMATFTE.continuationtype,'parameter')
    modelpar(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
if ~isempty(OCMATFTE.freeparameterindex)
    modelpar(OCMATFTE.freeparameterindex)=freepar(OCMATFTE.freeparametercoordinate);
end

arctime=OCMATFTE.initialarctimes;
%dtds=1;
%t=tcolmesh;
if strcmp(OCMATFTE.continuationtype,'time')
    arctime(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
diffarctime=diff(arctime);

t=diffarctime*tcolmesh+arctime(1);
dtds=diffarctime;

dxdt=dtds*OCMATFTE.canonicalsystem(t,depvar,modelpar);
dLdu=OCMATFTE.dlagrangefunctiondu(t,depvar,modelpar);
phi=OCMATFTE.complementaryslacknesscondition(t,depvar,modelpar);
rhs=[dxdt;dLdu;phi];
if OCMATFTE.objectivevaluecalc
    rhs(OCMATFTE.objectivevaluecoordinate,:)=dtds*OCMATFTE.objectivefunction(t,depvar,modelpar);
end

%-------------------------------------------------------------------------
function [J,Jpar]=daejac(tcolmesh,depvar,freepar,modelpar)
global OCMATFTE
arctime=OCMATFTE.initialarctimes;
if strcmp(OCMATFTE.continuationtype,'time')
    arctime(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
diffarctime=diff(arctime);

t=diffarctime*tcolmesh+arctime(1);
dtds=diffarctime;

J=OCMATFTE.dFDX;
Jpar=OCMATFTE.dFDPAR;
Jpartmp=[];
J(OCMATFTE.totalcoordinate,OCMATFTE.totalcoordinate)=[dtds*OCMATFTE.canonicalsystemjacobian(t,depvar,modelpar); ...
    OCMATFTE.dlagrangefunctiondujacobian(t,depvar,modelpar); ...
    OCMATFTE.complementaryslacknessconditionjacobian(t,depvar,modelpar)];
if OCMATFTE.objectivevaluecalc
    J(OCMATFTE.objectivevaluecoordinate,OCMATFTE.totalcoordinate)=dtds*OCMATFTE.objectivefunctionjacobian(t,depvar,modelpar);
end
if ~isempty(OCMATFTE.freeparameterindex)
    Jpartmp=[dtds*OCMATFTE.canonicalsystemparameterjacobian(t,depvar,modelpar); ...
        OCMATFTE.dlagrangefunctionduparameterjacobian(t,depvar,modelpar); ...
        OCMATFTE.complementaryslacknessconditionparameterjacobian(t,depvar,modelpar)];
    if OCMATFTE.objectivevaluecalc
        Jpartmp(OCMATFTE.objectivevaluecoordinate,:)=dtds*OCMATFTE.objectivefunctionparameterjacobian(t,depvar,modelpar);
    end
    Jpar(:,OCMATFTE.freeparametercoordinate)=Jpartmp(:,OCMATFTE.freeparameter);
end
if strcmp(OCMATFTE.continuationtype,'parameter')
    if isempty(Jpartmp)
        Jpartmp=[dtds*OCMATFTE.canonicalsystemparameterjacobian(t,depvar,modelpar); ...
            OCMATFTE.dlagrangefunctionduparameterjacobian(t,depvar,modelpar); ...
            OCMATFTE.complementaryslacknessconditionparameterjacobian(t,depvar,modelpar)];
        if OCMATFTE.objectivevaluecalc
            Jpartmp(OCMATFTE.objectivevaluecoordinate,:)=dtds*OCMATFTE.objectivefunctionparameterjacobian(t,depvar,modelpar);
        end
    end
    Jpar(:,OCMATFTE.continuationcoordinate)=Jpartmp(:,OCMATFTE.continuationindex);
end
if strcmp(OCMATFTE.continuationtype,'time')
    Jpartmp=Jpar(:,OCMATFTE.continuationcoordinate);
    Jpartmp([OCMATFTE.statecoordinate(:);OCMATFTE.costatecoordinate(:);])=OCMATFTE.canonicalsystem(t,depvar,modelpar);
    if OCMATFTE.objectivevaluecalc
        Jpartmp(OCMATFTE.objectivevaluecoordinate,:)=OCMATFTE.objectivefunction(t,depvar,modelpar);
    end
    Jpar(:,OCMATFTE.continuationcoordinate)=Jpartmp;
end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATFTE
if strcmp(OCMATFTE.continuationtype,'parameter')
    modelpar(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
if ~isempty(OCMATFTE.freeparameterindex)
    modelpar(OCMATFTE.freeparameterindex)=freepar(OCMATFTE.freeparametercoordinate);
end
arctime=OCMATFTE.initialarctimes;
if strcmp(OCMATFTE.continuationtype,'time')
    arctime(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end

timehorizon=arctime(end);

if strcmp(OCMATFTE.continuationtype,'initialstate')
    initstate=OCMATFTE.initialstate+freepar(OCMATFTE.continuationcoordinate)*OCMATFTE.targetvector;
else
    initstate=OCMATFTE.fixinitialstate;
end
resinit=OCMATFTE.bcinitial(arctime(1),depvara(:,1),OCMATFTE.fixinitialstatecoordinate,initstate,modelpar);
restrans=OCMATFTE.bctransversality(timehorizon,depvarb,modelpar);
if OCMATFTE.objectivevaluecalc
    OVal=OCMATFTE.salvagevalue(timehorizon,depvarb(:,end),modelpar);
    resinit=[resinit;depvara(OCMATFTE.objectivevaluecoordinate,1)-OVal];
end
res=[resinit;restrans(:)];
%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATFTE 
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tcolmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tcolmesh,coeff,tangent,varargin)
global OCMATCONT OCMATFTE

%----------------------------------------------------------------
function h=plotcontinuation(tcolmesh,coeff,tangent)
global OCMATFTE OCMATCONT
[tcolmesh,x,y,z1,freepar,modelpar]=drearr(tcolmesh,coeff,tangent);
arctime=OCMATFTE.initialarctimes;
if strcmp(OCMATFTE.continuationtype,'time')
    arctime(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
diffarctime=diff(arctime);

t=diffarctime*tcolmesh+arctime(1);
x=coeff2points(tcolmesh,coeff,'grid');
t=t(OCMATCONT.MESHDATA.tmeshidx);
h=OCMATFTE.plotcontinuation(t,x,modelpar);
%----------------------------------------------------------------
function idx=printcontinuation(tcolmesh,coeff,tangent)
global OCMATCONT OCMATFTE
idx=[];
if isempty(coeff)
    return
end
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.MESHDATA.continuationindex));
%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATFTE
varargout{1}=[];
%-------------------------------------------------------------
function [out, failed]=testfunc(id,tcolmesh,coeff,tangent)
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
        % Jacobian extended with bordering vectors v and w
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
function [out, failed]=targetvaluefunc(id,tcolmesh,coeff,tangent)
global OCMATCONT OCMATFTE

failed=[];
for ii=id
    switch ii
        case 1
            if strcmp(OCMATFTE.continuationtype,'initialstate')
                out=coeff(OCMATCONT.MESHDATA.continuationindex)-1;
            else
                out=coeff(OCMATCONT.MESHDATA.continuationindex)-OCMATFTE.targetvalue;
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function coeff=settargetvalue(id,tcolmesh,coeff,tangent)
global OCMATCONT OCMATFTE

for ii=id
    switch ii
        case 1
            if strcmp(OCMATFTE.continuationtype,'initialstate')
                coeff(OCMATCONT.MESHDATA.continuationindex)=1;
            else
                coeff(OCMATCONT.MESHDATA.continuationindex)=OCMATFTE.targetvalue;
            end
    end
end

%----------------------------------------------------------------
function out=formatsolution(tcolmesh,coeff,tangent)
global OCMATCONT OCMATFTE

[tcolmesh,x,y,z1,freepar,modelpar]=drearr(tcolmesh,coeff,tangent);
out=transform2nativematlab(tcolmesh,coeff,modelpar);
arctime=OCMATFTE.initialarctimes;
if strcmp(OCMATFTE.continuationtype,'time')
    arctime(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
out.arcinterval=arctime;

out.timehorizon=arctime(end);
out.modelparameter=modelpar;
out.modelname=OCMATFTE.modelname;

out.solver=func2str(OCMATCONT.newtonsolver);
out.solverinfo.method='dae';
out.solverinfo.collocationnumber=OCMATCONT.CollocationNumber;
out.solverinfo.collocationmethod=OCMATCONT.CollocationMethod;
out.solverinfo.daeorder=zeros(OCMATCONT.componentnumber,1);
out.solverinfo.daeorder(OCMATCONT.firstordercoordinate)=1;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tcolmesh=tcolmesh;
out.solverinfo.contclass='extremaldae4ft';
out.solverinfo.freeparameter=coeff(OCMATCONT.MESHDATA.freeparameterindex);
out.solverinfo.continuationtype=OCMATFTE.continuationtype;
out.solverinfo.fixinitialstatecoordinate=OCMATFTE.fixinitialstatecoordinate;
out.solverinfo.freeinitialstatecoordinate=OCMATFTE.freeinitialstatecoordinate;
out.solverinfo.fixendstatecoordinate=OCMATFTE.fixendstatecoordinate;
out.solverinfo.freeendstatecoordinate=OCMATFTE.freeendstatecoordinate;
out.solverinfo.objectivevaluecalc=OCMATFTE.objectivevaluecalc;
out.solverinfo.objectivevaluecoord=OCMATFTE.objectivevaluecoordinate;

%---------------------------------------------------------------------
function [failed,s]=process(id,tcolmesh,coeff,tangent,s)
global OCMATCONT OCMATFTE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.MESHDATA.continuationindex));
switch id
    case 1 % LP
end
failed=0;
%------------------------------------------------------------
function [S,L]=singmat

% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

S=[];
L=[];

%---------------------------------------------------------
function done
WorkspaceDone;
% ---------------------------------------------------------

function WorkspaceInit(tcolmesh,coeff,tangent)
% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------
function [tcolmesh,x,y,z1,freepar,modelpar]=drearr(tcolmesh,coeff,tangent)
global OCMATCONT OCMATFTE

modelpar=OCMATFTE.modelparameter;
freepar=coeff(OCMATCONT.MESHDATA.freeparameterindex); % parameter values
if strcmp(OCMATFTE.continuationtype,'parameter')
    modelpar(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
if ~isempty(OCMATFTE.freeparameterindex)
    modelpar(OCMATFTE.freeparameterindex)=freepar(OCMATFTE.freeparametercoordinate);
end
y=coeff(OCMATCONT.MESHDATA.ycoefficientidx); % y coefficients
z1=coeff(OCMATCONT.MESHDATA.firstorderzcoefficientidx); % z coefficients for the first order components
x=coeff2points_new(tcolmesh(1:OCMATCONT.CollocationNumber+1:end),coeff,'collocationgrid'); % returns the values of the states approximated by the collocation polynomials at the mesh with collocation points

% ------------------------------------------------------
function [tcolmesh,x,dxdt,freepar,modelpar]=residualrearr(tcolmesh,coeff,tangent,tres)
global OCMATCONT OCMATFTE

modelpar=OCMATFTE.modelparameter;
freepar=coeff(OCMATCONT.MESHDATA.freeparameterindex); % parameter values
if strcmp(OCMATFTE.continuationtype,'parameter')
    modelpar(OCMATFTE.continuationindex)=freepar(OCMATFTE.continuationcoordinate);
end
if ~isempty(OCMATFTE.freeparameterindex)
    modelpar(OCMATFTE.freeparameterindex)=freepar(OCMATFTE.freeparametercoordinate);
end
x=coeff2points(tcolmesh,coeff,'general',tres); % returns the values of the states approximated by the collocation polynomials at the mesh with collocation points
dxdt=coeff2derivatives(tcolmesh,coeff,'general',tres); % returns the values of the states approximated by the collocation polynomials at the mesh with collocation points

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATFTE
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATFTE=OCMATFTE;
try
    if contnum==1
        save([OCMATFTE.basicglobalvarfilename '4extremaldae4ft'],'MODELINFO')
    end
    save([OCMATFTE.basicresultfilename '4extremaldae4ft'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function failed=writelogfile(infostr)
global OCMATCONT OCMATFTE
persistent fid

failed=0;
if isempty(infostr)
    %fid=fopen(fullocmatfile(OCMATFTE.datapath(),['LogFile_' datestr(now,'mm_dd_yyyy_HH_MM_SS') '.log']),'w+');
    fclose('all');
    fid=fopen(fullocmatfile(OCMATFTE.datapath(),[OCMATFTE.modelname '_LogFile_' datestr(now,'mm_dd_yyyy') '.log']),'w+');
    if fid==-1
        failed=1;
        return
    end
    OCMATCONT.logfilehandle=fid;
elseif infostr==-1
    failed=fclose(fid);
else
    failed=fprintf(fid,infostr);
    %pause(0.1)
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATFTE

pathname=OCMATFTE.datapath();


%% functions for testing purposes
function testoperatoreq(tcolmesh,coeff,tangent,daefun,bcfun,icfun,res)
global OCMATCONT OCMATFTE

problemfile=['bvps2_' OCMATFTE.modelname '_probdef'];
ordnung=feval(problemfile,'orders');
m=OCMATCONT.CollocationNumber;
psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for i=1:m
        psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(i,OCMATCONT.CollocationPoint,ord);
    end
end
psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for i=1:m
        %evaluation of psi at the collocation points,0 and 1
        psival(ord,i,1:m)=polyval(psi(i,:,ord),OCMATCONT.CollocationPoint(1:m));
        psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
        psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
    end
end
tic
res0=functionFDF('F',problemfile,coeff,OCMATCONT.MESHDATA.tmesh,psival,OCMATCONT.FirstOrderCollocationPolynomial,OCMATCONT.CollocationPoint,[]);
t=toc;fprintf('elapsed time: %f\n',t)
max(abs(res-res0))


function Jbvpsuite=testfrechetderivative(tcolmesh,coeff,tangent,daefun,bcfun,icfun,daejacfun,bcjacfun,icfunjac,J,modelpar)
global OCMATCONT OCMATFTE MYTEST

% tic
numJacOpt.diffvar=2;
numJacOpt.vectvars=[];
Jnum=numjaccsd(@operatoreq,{tcolmesh,coeff,tangent,daefun,bcfun,icfun},numel(coeff),numJacOpt);
t=toc;fprintf('elapsed time: %f\n',t)
Jnum(end,:)=[];
Jnum(:,end)=[];
problemfile=['bvps2_' OCMATFTE.modelname '_probdef'];
ordnung=feval(problemfile,'orders');
m=OCMATCONT.CollocationNumber;
psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for i=1:m
        psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(i,OCMATCONT.CollocationPoint,ord);
    end
end
psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for i=1:m
        %evaluation of psi at the collocation points,0 and 1
        psival(ord,i,1:m)=polyval(psi(i,:,ord),OCMATCONT.CollocationPoint(1:m));
        psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
        psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
    end
end
MYTEST.par=modelpar;
tic
Jbvpsuite=functionFDF('DF',problemfile,coeff(1:end-1),tcolmesh(OCMATCONT.MESHDATA.tmeshidx),psival,OCMATCONT.FirstOrderCollocationPolynomial,OCMATCONT.CollocationPoint,[]);
t=toc;fprintf('elapsed time: %f\n',t)
J=J(1:size(Jbvpsuite,1),1:size(Jbvpsuite,2));
DJ=J-Jbvpsuite;
if max(abs(DJ(:)))>1e-10
    Jbvpsuite=full(Jbvpsuite);
    Jtmp=full(J);
    cols=1;
    rows=1;
    for ii=1:OCMATCONT.MESHDATA.meshNumber(1)
        if ii==1
            J_part=Jtmp(OCMATCONT.bcidx,:);
            Jnum_part=Jnum(OCMATCONT.bcidx,:);
            Jbvps_part=Jbvpsuite(OCMATCONT.bcidx,:);
        else
            siz=[OCMATCONT.MESHDATA.BasicJBlockSize(1)+OCMATCONT.MESHDATA.BasicCBlockSize(1) OCMATCONT.MESHDATA.BasicCBlockSize(2)];
            absoluteidx=rel2absidx(siz,[rows,cols],OCMATCONT.MESHDATA.continuationindex-1);
            J_part=reshape(Jtmp(absoluteidx),siz(1),[]);
            Jnum_part=reshape(Jnum(absoluteidx),siz(1),[]);
            Jbvps_part=reshape(Jbvpsuite(absoluteidx),siz(1),[]);
        end
        J_part
        Jnum_part
        Jbvps_part
        J_part-Jbvps_part
        answer=input('Interrupt the comparison (y)/n: ','s');
        while 1
            if isempty(answer)
                % default value 'y'
                answer='y';
            end
            if strcmpi(answer,'n')
                break
            elseif strcmpi(answer,'y')
                return
            end
        end
        if ii==1
            rows=rows+OCMATCONT.bcidx(end);
        else
            rows=rows+OCMATCONT.CollocationNumber*OCMATCONT.componentnumber;
            cols=cols+OCMATCONT.CollocationNumber*OCMATCONT.componentnumber;
        end
    end
end


function Psireturn=Psi(n,rho,nr)
% calculates the Lagrangepolynomials and its integrals
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;


