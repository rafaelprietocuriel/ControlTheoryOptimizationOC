function out=indifferencesolution4lcae()

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
if solutionindex==2
    arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate)];
else
    if OCMATINDIF.freeendtime
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.endtimecoordinate)];
    else
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.endtime];
    end
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
global OCMATINDIF OCMATCONT
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if solutionindex==2
    arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate)];
else
    if OCMATINDIF.freeendtime
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.endtimecoordinate)];
    else
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.endtime];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);

t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc);

J=OCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;

Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATINDIF.numarc(solutionindex)>1
    dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATINDIF.autonomous
        Jt=OCMATINDIF.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-transformedtimeshift-relarc+1)*Jt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-transformedtimeshift-relarc)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATINDIF.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-transformedtimeshift-relarc)*Jt);
        if OCMATINDIF.freeendtime && solutionindex==1
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.endtimecoordinate)=dxdt;
        else
            dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATINDIF.periodcoordinate)=dxdt;
        end
    end
else
    if OCMATINDIF.freeendtime && solutionindex==1
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATINDIF.endtimecoordinate)=dxdt;
    else
        dxdt=OCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATINDIF.periodcoordinate)=dxdt;
    end
end
Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATINDIF.freeparametercoordinate)=Jmodelpar(:,OCMATINDIF.freeparameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF
resinit=[];
resconnec=[];

depvaralc=depvara(:,OCMATINDIF.solutionindex==2);
depvarblc=depvarb(:,OCMATINDIF.solutionindex==2);

resper=OCMATINDIF.bclimitcycle(depvaralc,depvarblc);
resper=[resper; ...
    sum(OCMATINDIF.velocityvector.*(depvaralc(OCMATINDIF.velocitycoordinate,1)-OCMATINDIF.initialpoint))];


for ii=1:2
    actdepvara=depvara(:,OCMATINDIF.solutionindex==ii);
    actdepvarb=depvarb(:,OCMATINDIF.solutionindex==ii);
    switchtimes=freepar(OCMATINDIF.switchtimecoordinate{ii});
    arcarg=OCMATINDIF.arcarg{ii};
    %resinit=[resinit;depvara(end,OCMATINDIF.initcoord(ii))];
    for arc=1:OCMATINDIF.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
            OCMATINDIF.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
    end
    if ii==1
        resinit=depvaralc(OCMATINDIF.statecoordinate,1)-actdepvara(OCMATINDIF.statecoordinate,1);

        hatx=freepar(OCMATINDIF.equilibriumcoordinate);
        resequilibrium=OCMATINDIF.canonicalsystem(0,hatx,modelpar,arcarg(end));
        Jac=OCMATINDIF.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
        Y=freepar(OCMATINDIF.Ycoordinate);
        asymptoticmatrix=OCMATINDIF.Q0*[-Y';OCMATINDIF.Id];
        resricatti=ricatti(Y,Jac,OCMATINDIF.Q0,OCMATINDIF.subspacedim);
        resasym=OCMATINDIF.bcasymptotic(actdepvarb(:,end),asymptoticmatrix,hatx);
        if OCMATINDIF.freeendtime
            resasym=[resasym; ...
                sqrt(sum((hatx-actdepvarb(:,end)).^2))-OCMATINDIF.distance(ctr)];
        end
        residpt=OCMATINDIF.bcindifference([actdepvara(:,1) depvaralc(:,1)],modelpar,[OCMATINDIF.arcarg{1}(1) OCMATINDIF.arcarg{2}(1)],[1 2]);
        if length(OCMATINDIF.freeparameterindex)==1
            residpt=residpt-freepar(end)*OCMATINDIF.hamiltoniandifference;
        end
    end
end
res=[resinit;resasym;resconnec;residpt;resricatti;resper;resequilibrium];
%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
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
    if solutionindex==2
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate)];
    else
        if OCMATINDIF.freeendtime
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.endtimecoordinate)];
        else
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.endtime];
        end
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
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATINDIF.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
%global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
fprintf(1,' Continuation parameter: %g\n',coeff(end));

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
            if length(OCMATINDIF.freeparameterindex)==1
                out=coeff(OCMATCONT.HE.contparametercoord);

            else
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
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

out=transform2nativematlab(tmesh,coeff,OCMATINDIF.parametervalue);
out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
    OCMATCONT.HE.TIMEDDATA.rightarcindex];
arctime=[];
for solutionindex=1:2
    if solutionindex==2
        arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.periodcoordinate)];
    else
        if OCMATINDIF.freeendtime
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' freepar(OCMATINDIF.endtimecoordinate)];
        else
            arctime=[0 freepar(OCMATINDIF.switchtimecoordinate{solutionindex}).' OCMATINDIF.endtime];
        end
    end
    out.solverinfo.arcinterval{solutionindex}=arctime;
    out.solverinfo.timehorizon(solutionindex)=arctime(end);
end
out.arcinterval=arctime;
out.arcarg=OCMATCONT.HE.arcarg;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution4lcae';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATINDIF.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoordinate=OCMATINDIF.switchtimecoordinate;
out.solverinfo.periodcoordinate=OCMATINDIF.periodcoordinate;
if OCMATINDIF.freeendtime
    out.solverinfo.endtimecoordinate=OCMATINDIF.endtimecoordinate;
end
out.solverinfo.Ycoordinate=OCMATINDIF.Ycoordinate;
out.solverinfo.equilibriumcoordinate=OCMATINDIF.equilibriumcoordinate;
out.solverinfo.freeparametercoordinate=OCMATINDIF.freeparametercoordinate;
out.solverinfo.freeparameterindex=OCMATINDIF.freeparameterindex;
out.solverinfo.solutionindex=OCMATINDIF.solutionindex;


%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
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


%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

% ------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF

modelpar=OCMATINDIF.parametervalue;

y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
z=[];
freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
modelpar(OCMATINDIF.freeparameterindex)=freepar(OCMATINDIF.freeparametercoordinate);
%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATINDIF OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATINDIF=OCMATINDIF;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4lcae'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4lcae'],'sout','bvpout')
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
    [R11,R12,R21,R22]=RicattiCoeff(Q0,J,dim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
