function out=extremalp2ep()

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
out{12}=@objectivederivative;
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

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});


%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
modelpar(OCMATAE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATAE.objectivevaluecalc
    dxdt(OCMATAE.objectivevaluecoord,:)=dtds*OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
end
if OCMATAE.exogenousfunction
    dxdt(OCMATAE.exogenousdynamicscoordinate,:)=dtds*OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
modelpar(OCMATAE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
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
if OCMATAE.exogenousfunction
    J=[[J zeros(size(J,1),1)];dtds*OCMATAE.exogenousjacobian(t,depvar,modelpar,arcarg)];
end

Jpar=zeros(OCMATAE.ODEcoord(end),OCMATCONT.HE.numparameter);
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
    if OCMATAE.exogenousfunction
        dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
    end
    if arc<OCMATCONT.HE.numarc
        Jpar(OCMATAE.ODEcoord,OCMATAE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        if arc>1
            Jpar(OCMATAE.ODEcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
    else
        Jpar(OCMATAE.ODEcoord,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        if OCMATAE.movinghorizon
            %dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATAE.movinghorizoncoord)=dxdt;
        end
    end
else
    if OCMATAE.movinghorizon
        dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
        if OCMATAE.exogenousfunction
            dxdt(OCMATAE.exogenousdynamicscoordinate,:)=OCMATAE.exogenousdynamics(t,depvar,modelpar,arcarg);
        end
        Jpar(:,OCMATAE.movinghorizoncoord)=dxdt;
    end
end
Jmodelpar=dtds*OCMATAE.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
if OCMATAE.exogenousfunction
    Jexogenousmodelpar=OCMATAE.exogenousparameterjacobian(t,depvar,modelpar,arcarg);
    Jmodelpar=[Jmodelpar;Jexogenousmodelpar];
end
Jpar(:,OCMATCONT.HE.numparameter)=Jmodelpar(:,OCMATAE.varyparameterindex);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
modelpar(OCMATAE.varyparameterindex)=freepar(OCMATCONT.HE.numparameter);
switchtimes=freepar(OCMATAE.switchtimecoord);

resconnec=[];
resuser=[];
hatx=freepar(OCMATCONT.HE.equilibriumcoord);
if ~isempty(OCMATAE.fixcoordinate)
    hatx(OCMATAE.varcoordinate)=hatx;
    hatx(OCMATAE.fixcoordinate)=OCMATAE.fixcoordinatevalue;
end
Y=freepar(OCMATCONT.HE.Ycoord);
asymptoticmatrix=OCMATAE.Q0*[-Y';OCMATAE.Id];
Jac=OCMATAE.canonicalsystemjacobian(0,hatx,modelpar,OCMATAE.arcid4ep);
if ~isempty(OCMATAE.excludecoordinate4ep)
    Jac(:,OCMATAE.excludecoordinate4ep)=[];
    Jac(OCMATAE.excludecoordinate4ep,:)=[];
end
resequilibrium=OCMATAE.equilibrium(hatx,modelpar,OCMATAE.arcid4ep);
%resequilibrium=OCMATAE.canonicalsystem(0,hatx,modelpar,OCMATCONT.HE.arcarg(end));
resequilibrium(OCMATAE.fixcoordinate)=[];
resricatti=ricatti(Y,Jac);
if ~OCMATAE.followequilibrium
    resinit=OCMATAE.bcinitial(depvara,OCMATAE.initialcoordinate,OCMATAE.initialstate,modelpar,OCMATCONT.HE.arcarg(1));
else
    resinit=OCMATAE.bcinitial(depvara,OCMATAE.initialcoordinate,hatx(OCMATAE.initialcoordinate),modelpar,OCMATCONT.HE.arcarg(1));
end
if OCMATAE.exogenousfunction
    resinit=[resinit; ...
        depvara(OCMATAE.exogenousdynamicscoordinate,1)-OCMATAE.exogenousinitialstates];
end

if OCMATAE.objectivevaluecalc
    resinit=[resinit;depvara(end,1)];
end
if ~isempty(OCMATAE.userbc)
    if ~OCMATAE.movinghorizon
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
    end
    resuser=OCMATAE.userbc(arctime(end),depvar,depvarb,par,OCMATCONT.HE.arcarg(1));
end
try
    resasym=OCMATAE.bcasymptotic(depvarb,asymptoticmatrix,hatx);
catch
    resasym=OCMATAE.bcasymptotic(depvarb,asymptoticmatrix,hatx,modelpar,OCMATCONT.HE.arcarg(end),depvara);
end
if OCMATAE.movinghorizon
    yend=depvarb(:,end);
    yend=yend(1:length(hatx));
    if ~isempty(OCMATAE.excludecoordinate4ep)
        yend(OCMATAE.excludecoordinate4ep,:)=[];
    end
    resasym=[resasym; ...
        sqrt(sum((yend-hatx).^2))-OCMATAE.distance];
end

for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
end

res=[resinit;resequilibrium;resricatti;resconnec;resasym;resuser];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvaral,depvarar,depvarbl,depvarbr,modelpar,freepar,contval,arc)
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

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    if ~OCMATAE.movinghorizon
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
    else
        arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
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
    if strcmp(OCMATAE.pathtype,'u')
        violationmat=-diffarctime(arc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    else
        violationmat=diffarctime(arc)<-OCMATCONT.OPTIONS.zerotimedifftolerance;
    end
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
if ~OCMATAE.movinghorizon
    hatx=freepar(OCMATCONT.HE.equilibriumcoord);
    yend=y(:,end);
    if ~isempty(OCMATAE.excludecoordinate4ep)
        yend(OCMATAE.excludecoordinate4ep,:)=[];
    end
    
    violationmat=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx)<0;
    if violationmat
        counter=counter+1;
        cols=size(y,2);
        infoS(counter).arcarg=arcarg;
        infoS(counter).arcnum=arc;
        infoS(counter).rows='maxdistance';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=norm(norm(yend-hatx));
        infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(yend-hatx);
        b=min([b infoS(counter).minval]);
    end
end

%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
if ~OCMATAE.movinghorizon
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
end
diffarctime=diff(arctime);
for arc=1:numel(OCMATCONT.HE.arcarg)
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(arc)-diffarctime(arc)*(arc-1));
end
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATAE.plotcontinuation(t,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT OCMATAE
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
fprintf(1,' Continuation parameter: %2.9g\n',coeff(OCMATCONT.HE.contparametercoord));
if OCMATAE.findoptimalparameter
    if OCMATAE.objectivevaluecalc
        tpar=tangent(end);
        tangent=tangent/tpar;
        tangent=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
        out=tangent(OCMATAE.objectivevaluecoord,end);
    else
        out=objectivederivative(tmesh,coeff,tangent);
    end
    fprintf(1,' Derivative value      : %g\n',out);
elseif ~isempty(OCMATAE.targetfunction)
    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
    out=OCMATAE.targetfunction(t(1),y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex),y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex),modelpar,OCMATCONT.HE.arcarg);
    fprintf(1,' Difference to target value      : %g\n',out);

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
            if OCMATAE.findoptimalparameter
                if OCMATAE.objectivevaluecalc
                    tpar=tangent(end);
                    tangent=tangent/tpar;
                    tangent=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
                    out=tangent(OCMATAE.objectivevaluecoord,end);
                else
                    out=objectivederivative(tmesh,coeff,tangent);
                end
                OCMATAE.derivative=out;
            elseif ~isempty(OCMATAE.targetfunction)
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
                out=OCMATAE.targetfunction(t(1),y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex),y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex),modelpar,OCMATCONT.HE.arcarg);
            else
                out=OCMATAE.targetparametervalue-coeff(OCMATCONT.HE.contparametercoord);
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATAE OCBVP

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=[1;length(out.x)];
end
if ~OCMATAE.movinghorizon
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
else
    out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATAE.movinghorizoncoord)];
    out.solverinfo.distance=OCMATAE.distance;
    out.solverinfo.movinghorizoncoord=OCMATAE.movinghorizoncoord;
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
out.timehorizon=OCMATAE.truncationtime;
out.modelparameter=modelpar;
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
out.solverinfo.subspacedim=OCMATAE.subspacedim;
out.solverinfo.orthspacedim=OCMATAE.orthspacedim;
out.solverinfo.qbasis=OCMATAE.Q0;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;
out.solverinfo.inftimetransformation=OCMATAE.inftimetransformation;
out.solverinfo.switchtimecoord=OCMATAE.switchtimecoord;
if OCMATAE.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATAE.objectivevaluecoord;
end
out.solverinfo.findoptimalparameter=OCMATAE.findoptimalparameter;
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.solverinfo.yp=out.yp;
        out.solverinfo.ypmid=out.ypmid;
        out=rmfield(out,{'yp','ypmid'});
    case 'bvp4c'
        out.solverinfo.yp=out.yp;
        out=rmfield(out,'yp');
end
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

modelpar=OCMATAE.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp4c'}
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
[U,S,V]=svd(OCMATAE.Q0(:,1:OCMATAE.subspacedim)+OCMATAE.Q0(:,OCMATAE.subspacedim+1:end)*Y);
OCMATAE.Q0=U;
OCMATAE.Y=zeros(OCMATAE.orthspacedim,OCMATAE.subspacedim);

freepar(OCMATCONT.HE.Ycoord)=OCMATAE.Y;

switch OCMATCONT.bvpmethod
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
function out=ricatti(Y,J)
global OCMATAE
out=[];

if  ~isempty(Y)
    % Riccati blocks from unstable eigenspace
    [R11,R12,R21,R22]=RicattiCoeff(OCMATAE.Q0,J,OCMATAE.subspacedim);
    out=R22*Y-Y*R11+R21-Y*R12*Y;
end
out=out(:);

function out=objectivederivative(tmesh,coeff,tangent)
global OCMATCONT OCMATAE

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
tangenty=tangent(OCMATCONT.HE.DDATA.meshvalcoord);
tangentp=tangent(OCMATCONT.HE.contparametercoord); % parameter values


arcarg=OCMATCONT.HE.arcarg(1);
tangent1=[tangenty(:,1);tangentp];
tangent1=tangent1/norm(tangent1);
[dHdx,dHdp]=OCMATAE.hamiltonianderivative(t(1),y(:,1),modelpar,arcarg);
out=dHdx(:).'*tangent1(1:end-1)+dHdp(OCMATAE.varyparameterindex)*tangent1(end);

