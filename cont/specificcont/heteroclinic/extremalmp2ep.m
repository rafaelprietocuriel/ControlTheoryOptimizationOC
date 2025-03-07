function out=extremalmp2ep()

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
global OCMATAE OCMATCONT OCBVP
solutionindex=OCMATAE.solutionindex(arc);
modelpar(OCMATAE.parameterindex)=freepar(OCMATAE.parametervaluecoord);

if solutionindex>1
    relarc=arc-OCMATAE.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if strcmp(OCMATAE.limitsettype,'e')
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' OCMATAE.truncationtime(solutionindex)];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' freepar(OCMATAE.truncationtimecoord{solutionindex})];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATAE.numarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
solutionindex=OCMATAE.solutionindex(arc);

if solutionindex>1
    relarc=arc-OCMATAE.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if strcmp(OCMATAE.limitsettype,'e')
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' OCMATAE.truncationtime(solutionindex)];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' freepar(OCMATAE.truncationtimecoord{solutionindex})];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
if solutionindex>1
    transformedtimeshift=OCMATAE.numarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATAE.numarc(solutionindex)>1
    dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATAE.autonomous
        Jt=OCMATAE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if relarc<OCMATAE.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-relarc)*Jt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
    end
end
Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATAE.parametervaluecoord)=Jmodelpar(:,OCMATAE.parameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
resinit=[];
resasym=[];
resconnec=[];
resequilibrium=[];
resricatti=[];

modelpar(OCMATAE.parameterindex)=freepar(OCMATAE.parametervaluecoord);

arcarg=OCMATCONT.HE.arcarg(OCMATAE.arccoord{1});
hatx=freepar(OCMATCONT.HE.equilibriumcoord);
resequilibrium=[resequilibrium;OCMATAE.canonicalsystem(0,hatx,modelpar,arcarg(end))];
for ii=1:OCMATAE.multorder
    switchtimes=freepar(OCMATAE.switchtimecoord{ii});
    arcarg=OCMATCONT.HE.arcarg(OCMATAE.arccoord{ii});
    resinit=[resinit; ...
        OCMATAE.bcinitial(depvara(:,OCMATAE.initcoord(ii)),OCMATAE.initialcoordinate,OCMATAE.initialstate,modelpar,OCMATCONT.HE.arcarg(1));];
    if OCMATAE.simple
        asymptoticmatrix=OCMATAE.asymptoticmatrix{ii};
    else
        Y=freepar(OCMATCONT.HE.Ycoord{ii});
        asymptoticmatrix=OCMATAE.Q0{ii}*[-Y';OCMATAE.Id{ii}];
        Jac=OCMATAE.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
        resricatti=[resricatti;ricatti(Y,Jac,OCMATAE.Q0{ii},OCMATAE.subspacedim{ii})];
    end
    resasym=[resasym;OCMATAE.bcasymptotic(depvarb(:,OCMATAE.cumsumnumarc(ii)),asymptoticmatrix,hatx)];
    for arc=1:OCMATAE.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATAE.reset(depvara(:,OCMATAE.arccoord{ii}),depvarb(:,OCMATAE.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATAE.edge{ii},arc); ...
            OCMATAE.guard(depvara(:,OCMATAE.arccoord{ii}),depvarb(:,OCMATAE.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATAE.edge{ii},arc)];
    end
end
res=[resinit;resasym;resconnec;resequilibrium;resricatti];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATAE.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATAE.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATAE.numarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' OCMATAE.truncationtime(solutionindex)];
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
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
    if OCMATAE.stableflag{solutionindex}
        violationmat=diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    else
        violationmat=-diffarctime(relarc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    end
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
global OCMATAE OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
if length(sol.x)~=size(sol.y,2)
    sol
end
% clear possible persistent variable
h=OCMATAE.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATAE OCBVP
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
fprintf(1,' Continuation parameter: %g\n',coeff(end));
fprintf(1,' Difference initial state: %g\n',norm(y(:,OCBVP.Lidx(OCMATAE.initcoord(1)))-y(:,OCBVP.Lidx(OCMATAE.initcoord(2)))));

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
for ii=id
    switch ii
        case 1
            out=OCMATAE.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
            %out=tangent(end);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATAE OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[];
for solutionindex=1:2
    %solutionindex=OCMATAE.solutionindex(arc);
    if strcmp(OCMATAE.limitsettype,'e')
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' OCMATAE.truncationtime(solutionindex)];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord{solutionindex}).' freepar(OCMATAE.truncationtimecoord{solutionindex})];
    end
    out.arcinterval=[out.arcinterval arctime];
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
out.timehorizon=OCMATAE.truncationtime;
%out.modelparameter=OCMATAE.parametervalue;
%out.modelparameter(OCMATAE.parameterindex)=freepar(OCMATAE.parametervaluecoord);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='heteroclinic';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATAE.inftimetransformation;
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATAE.switchtimecoord;
out.solverinfo.truncationtimecoord=OCMATAE.truncationtimecoord;
out.solverinfo.equilibriumcoord=OCMATCONT.HE.equilibriumcoord;
out.solverinfo.Ycoord=OCMATCONT.HE.Ycoord;
out.solverinfo.subspacedim=OCMATAE.subspacedim;
out.solverinfo.orthspacedim=OCMATAE.orthspacedim;
out.solverinfo.qbasis=OCMATAE.Q0;
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
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        % reduce Jacobian to size without the continuation parameter
        J(:,OCMATCONT.HE.contparametercoord)=[];
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        q=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        q=q/norm(q);
        p=Q(:,end);
        p=p/norm(p);
        s.data.phi=q(:);
        s.data.psi=p(:);
        s.data.DFDX=J;
        s.data.laecoefficient=calchnf_LP(tmesh,coeff,tangent,@operatoreq,q,p,J);
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
global OCMATCONT OCMATAE

%numarc=OCMATCONT.HE.numarc;
%domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATAE.parametervalue;
if length(OCMATAE.parameterindex)==2
    modelpar(OCMATAE.parameterindex)=coeff(OCMATCONT.HE.contparametercoord+(-1:0));
else
    modelpar(OCMATAE.parameterindex)=coeff(OCMATCONT.HE.contparametercoord-1);
end
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
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
        save([OCMATAE.basicglobalvarfilename '4extremalmp2ep'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremalmp2ep'],'sout','bvpout')
catch
    failed=1;
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
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATAE

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

for ii=1:OCMATAE.multorder
    Y=freepar(OCMATCONT.HE.Ycoord{ii});
    OCMATCONT.adapted = 1;
    %
    [U,S,V]=svd(OCMATAE.Q0{ii}(:,1:OCMATAE.subspacedim{ii})+OCMATAE.Q0{ii}(:,OCMATAE.subspacedim{ii}+1:end)*Y);
    OCMATAE.Q0{ii}= U;
    OCMATAE.Y{ii}=zeros(OCMATAE.orthspacedim{ii},OCMATAE.subspacedim{ii});

    freepar(OCMATCONT.HE.Ycoord{ii})=OCMATAE.Y{ii};
    switch OCMATCONT.bvpmethod
        case {'bvp6c','bvp4c'}
            coeff=[y(:);freepar];
        otherwise
    end
end
flag = 1;

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
