function out=extremalp2epuser()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@ode;
out{4}{2}=@bc;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;
out{5}{3}=@bctransversality;
out{5}{4}=@equilibrium;
out{5}{5}=@ricatti;

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
if ~OCBVP.multiarccalc
    dxdt=zeros(size(depvar));
    OCBVP.multiarccalc=1;
    for ii=1:OCMATCONT.HE.numarc
        depvarint=depvar(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:);
        dxdt(OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii),:)=ode(s,depvarint,ii,freepar,modelpar);
    end
    OCBVP.multiarccalc=0;
    return
end
modelpar(OCMATAE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if OCMATAE.inftimetransformation && arc==OCMATCONT.HE.numarc
    t=arctime(arc)+log(arc-s)/OCMATAE.inftimetransformation;
    dtds=-1./(OCMATAE.inftimetransformation*(arc-s));
else
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1)); 
    dtds=diffarctime(arc);
end
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    J=zeros(OCBVP.n);
    Jpar=zeros(OCBVP.n,OCBVP.npar);
    OCBVP.multiarccalc=1;
    for ii=1:OCMATCONT.HE.numarc
        idx=OCMATCONT.HE.arcrowindex(1,ii):OCMATCONT.HE.arcrowindex(2,ii);
        depvarint=depvar(idx,:);
        [J(idx,idx) Jpar(idx,1:OCBVP.npar)]=odejac(s,depvarint,ii,freepar,modelpar);
    end
    OCBVP.multiarccalc=0;
    return
end
arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
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
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATCONT.HE.numarc>1
    dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord(arc))=dxdt;
        if arc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord(arc-1))=-dxdt;
        end
    else
        if OCMATAE.inftimetransformation
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=0;
            %Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,arc-1)=-1./((exp(OCMATAE.inftimetransformation*arctime(arc))*(arc-s)))*dxdt;
        else
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATAE.switchtimecoord(arc-1))=-dxdt;
        end
    end
end
Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATAE.varyparameterindex);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    res=bc(depvara,depvarb,freepar,modelpar);
    OCBVP.multiarccalc=0;
    return
end
modelpar(OCMATAE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
switchtimes=freepar(OCMATAE.switchtimecoord);
if OCMATAE.stateconstraint
    jump=zeros(1,OCMATCONT.HE.numarc);
    jump(OCMATAE.jumpcostateindex)=freepar(OCMATAE.jumpcostatecoord);
end

resconnec=[];
hatx=freepar(OCMATCONT.HE.equilibriumcoord);
Y=freepar(OCMATCONT.HE.Ycoord);
asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
Jac=OCMATAE.canonicalsystemjacobian(0,hatx,modelpar,OCMATCONT.HE.arcarg(end));
resequilibrium=OCMATAE.canonicalsystem(0,hatx,modelpar,OCMATCONT.HE.arcarg(end));
resricatti=ricatti(hatx,Y,modelpar,OCMATCONT.HE.arcarg,OCMATAE.stableflag,OCMATAE.Q0,Jac);
%resinit=OCMATAE.bcinitial(depvara,OCMATAE.initialcoordinate,OCMATAE.initialstate,modelpar,OCMATCONT.HE.arcarg(1));
resinit=OCMATAE.bcinitial(depvara,OCMATAE.initialcoordinate,hatx,modelpar,OCMATCONT.HE.arcarg(1));
resasym=OCMATAE.bcasymptotic(depvarb,asymptoticmatrix,hatx);
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
res=[resinit;resequilibrium;resricatti;resconnec;resasym];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,contval,arc)
global OCMATAE OCMATCONT OCBVP
if ~OCBVP.multiarccalc
    OCBVP.multiarccalc=1;
    depvara=depvara(OCMATCONT.HE.bcindex);
    depvarb=depvarb(OCMATCONT.HE.bcindex);
    [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar);
    Ja(:,OCBVP.nBCs-OCBVP.nparmc+1:OCBVP.nBCs)=[];
    Jb(:,OCBVP.nBCs-OCBVP.nparmc+1:OCBVP.nBCs)=[];
    OCBVP.multiarccalc=0;
    return
end
Jalinit=[];
Jarinit=[];
Jblinit=[];
Jbrinit=[];
Jparinit=[];
Jalasym=[];
Jarasym=[];
Jblasym=[];
Jbrasym=[];
Jparasym=[];
Jalequilib=[];
Jarequilib=[];
Jblequilib=[];
Jbrequilib=[];
Jparequilib=[];
Jalricatti=[];
Jarricatti=[];
Jblricatti=[];
Jbrricatti=[];
Jparricatti=[];
Jalguard=[];
Jarguard=[];
Jblguard=[];
Jbrguard=[];
Jparguard=[];
Jalreset=[];
Jarreset=[];
Jblreset=[];
Jbrreset=[];
Jparreset=[];

arcarg=OCMATCONT.HE.arcarg(arc);
if isempty(OCMATCONT.HE.edge)
    Y=freepar(OCMATCONT.HE.Ycoord);
    asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
    hatx=freepar(OCMATCONT.HE.equilibriumcoord);
    v=-OCMATAE.Q0'*(depvarbl-hatx);
    v=v(OCMATAE.TemplateBCJacAsymptoticYIndex);
    Jac=OCMATAE.canonicalsystemjacobian(0,hatx,modelpar,arcarg);
    Hess=OCMATAE.canonicalsystemhessian(0,hatx,modelpar,arcarg);
    Jacpar=OCMATAE.canonicalsystemparameterjacobian(0,hatx,modelpar,arcarg);
    Hesspar=OCMATAE.canonicalsystemparameterhessian(0,hatx,modelpar,arcarg);
    Hesspar=Hesspar(:,:,OCMATAE.varyparameterindex);
    % the order of the variables is 
    % [y0 y1 par] with par=[hatx Y contpar]
    Jarinit=OCMATAE.BCJacinitialY0;
    Jblinit=OCMATAE.BCJacinitialYN;
    Jparinit=OCMATAE.BCJacinitialP;
    
    Jarasym=OCMATAE.TemplateBCJacAsymptoticY0;
    Jblasym=asymptoticmatrix';
    Jparasym=OCMATAE.TemplateBCJacAsymptoticP;
    Jparasym(:,OCMATCONT.HE.equilibriumcoord)=-asymptoticmatrix';
    JY=OCMATAE.TemplateBCJacAsymptoticY;
    JY(OCMATAE.TemplateBCJacAsymptoticNoneZeroIndex)=v(OCMATAE.TemplateBCJacAsymptoticvNoneZeroIndex);
    Jparasym(:,OCMATCONT.HE.Ycoord)=JY;
    
    Jarequilib=OCMATAE.TemplateBCJacEquilibriumY0;
    Jblequilib=OCMATAE.TemplateBCJacEquilibriumYN;
    Jparequilib=OCMATAE.TemplateBCJacEquilibriumP;
    Jparequilib(:,OCMATCONT.HE.equilibriumcoord)=Jac;
    Jparequilib(:,end)=Jacpar(:,OCMATAE.varyparameterindex);
    
    Jarricatti=OCMATAE.TemplateBCJacRicattiY0;
    Jblricatti=OCMATAE.TemplateBCJacRicattiYN;
    Jparricatti=OCMATAE.TemplateBCJacRicattiP;
    
    % Ricatti with respect to hatx
    for ii=1:numel(OCMATCONT.HE.equilibriumcoord)
        Jparricatti(:,OCMATCONT.HE.equilibriumcoord(ii))=ricatti(hatx,Y,modelpar,arcarg,OCMATAE.stableflag,OCMATAE.Q0,Hess(:,:,ii));
    end
    [R11,R12,R21,R22]=RicattiCoeff(OCMATAE.Q0,Jac,OCMATAE.dimSubSpace);
    R11=R11(:,:,ones(1,OCMATAE.numY));
    D1=OCMATAE.TemplateBCJacRicattiC;
    D1(OCMATAE.TemplateBCJacRicattiC1NoneZeroIndex)=-R11(OCMATAE.TemplateBCJacRicattiB1NoneZeroIndex);
    R22=R22(:,:,ones(1,OCMATAE.numY));
    D2=OCMATAE.TemplateBCJacRicattiC;
    D2(OCMATAE.TemplateBCJacRicattiC2NoneZeroIndex)=R22(OCMATAE.TemplateBCJacRicattiB2NoneZeroIndex);
    B=R12*Y;
    B=B(:,:,ones(1,OCMATAE.numY));
    D3=OCMATAE.TemplateBCJacRicattiC;
    D3(OCMATAE.TemplateBCJacRicattiC1NoneZeroIndex)=-B(OCMATAE.TemplateBCJacRicattiB1NoneZeroIndex);
    B=Y*R12;
    B=B(:,:,ones(1,OCMATAE.numY));
    D4=OCMATAE.TemplateBCJacRicattiC;
    D4(OCMATAE.TemplateBCJacRicattiC2NoneZeroIndex)=-B(OCMATAE.TemplateBCJacRicattiB2NoneZeroIndex);
    Jparricatti(:,OCMATCONT.HE.Ycoord(:))=D1+D2+D3+D4;
    Jparricatti(:,OCMATCONT.HE.numparameter)=ricatti(hatx,Y,modelpar,arcarg,OCMATAE.stableflag,OCMATAE.Q0,Hesspar);
elseif arc==1
    [Jalinit Jarinit Jblinit Jbrinit Jparinit]=OCMATAE.bcjacobianinitial(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,arcarg,OCMATCONT.HE.edge(arc,:),OCMATAE.initialcoordinate,OCMATAE.initialstate);
    [Jalguard Jarguard Jblguard Jbrguard Jparguard]=OCMATAE.jacobianguard(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,arcarg,OCMATCONT.HE.edge(arc,:));
elseif arc==OCMATCONT.HE.numarc
    [Jalreset Jarreset Jblreset Jbrreset Jparreset]=OCMATAE.jacobianreset(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,arcarg,OCMATCONT.HE.edge(arc-1,:));
    [Jalasym Jarasym Jblasym Jbrasym Jparasym]=OCMATAE.bcjacobianasymptotic(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,arcarg,OCMATCONT.HE.edge(arc-1,:),OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint);
else
    [Jalreset Jarreset Jblreset Jbrreset Jparreset]=OCMATAE.jacobianreset(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,arcarg,OCMATCONT.HE.edge(arc-1,:));
    [Jalguard Jarguard Jblguard Jbrguard Jparguard]=OCMATAE.jacobianguard(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,arcarg,OCMATCONT.HE.edge(arc,:));
end
Jal=[Jalinit;Jalreset;Jalguard;Jalasym;Jalequilib;Jalricatti];
Jar=[Jarinit;Jarreset;Jarguard;Jarasym;Jarequilib;Jarricatti];
Jbl=[Jblinit;Jblreset;Jblguard;Jblasym;Jblequilib;Jblricatti];
Jbr=[Jbrinit;Jbrreset;Jbrguard;Jbrasym;Jbrequilib;Jbrricatti];
Jpar=[Jparinit;Jparreset;Jparguard;Jparasym;Jparequilib;Jparricatti];

%----------------------------------------------------------------
function [res,varargout]=residual(tmesh,coeff,rhs,odefun)
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        [res resindex]=residual_sbvpoc(t,y,z,freepar,modelpar,@ode,rhs);
        varargout{1}=resindex;
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
    case 'bvp4c'
        res=residual_bvp4c(t,y,freepar,modelpar,rhs(1:OCMATCONT.HE.numdvariablesmc),odefun);
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
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
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
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATAE.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));

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
            out=OCMATAE.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATAE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=[1;length(out.x)];
end
out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
if OCMATAE.inftimetransformation
    out.timehorizon=inf;
else
    out.timehorizon=OCMATAE.truncationtime;
end
out.modelparameter=OCMATAE.parametervalue;
out.modelparameter(OCMATAE.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tmesh=tmesh;
out.solverinfo.tangent=tangent;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremalp2ep';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATAE.inftimetransformation;
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;

out.solverinfo.equilibriumcoord=OCMATCONT.HE.equilibriumcoord;
out.solverinfo.Ycoord=OCMATCONT.HE.Ycoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.solverinfo.inftimetransformation=OCMATAE.inftimetransformation;
out.solverinfo.switchtimecoord=OCMATAE.switchtimecoord;
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
if OCMATAE.stateconstraint
    out.solverinfo.jumpcostatecoord=OCMATAE.jumpcostatecoord;
end
%OCMATAE.initialstate=deval(tmesh,coeff,tangent,10/OCMATAE.truncationtime);
% add solver method specific information to make it consistent with MATLAB
% syntax
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
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

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATAE.parametervalue;
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

modelpar(OCMATAE.varyparameterindex)=coeff(OCMATCONT.HE.contparametercoord);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATAE

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff); 
Y=freepar(OCMATCONT.HE.Ycoord);
OCMATCONT.adapted = 1;
% 
Q=OCMATAE.Q0;
[U,S,V]=svd(Q*[eye(size(Y,1));Y']);
OCMATAE.Q0= U;
OCMATAE.Y=zeros(OCMATAE.numunstable,OCMATAE.numstable);
freepar(OCMATCONT.HE.Ycoord)=OCMATAE.Y;
switch OCMATCONT.bvpmethod
    case 'bvp5c'
        Y=[];
        for arc=1:OCMATCONT.HE.numarc
            Y=[Y y(OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc),:)];
        end
        Y=[Y;freepar(:,ones(1,size(Y,2)))];
        coeff=Y(:);
    case {'bvp6c','bvp4c'}
        coeff=[y(:);freepar];
    otherwise
end
flag = 1;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATAE OCBVP
failed=0; 
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATAE=OCMATAE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATAE.basicglobalvarfilename '4extremalp2ep'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremalp2ep'],'sout','bvpout')
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

%-----------------------------------------------------------------
function out=ricatti(xhat,Y,modelpar,arcarg,stableflag,Q,J)
global OCMATAE
out=[];

if stableflag
    dimInvSubSpace=OCMATAE.numstable;
else
    dimInvSubSpace=OCMATAE.numunstable;
end
if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(Q,J,dimInvSubSpace);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);
