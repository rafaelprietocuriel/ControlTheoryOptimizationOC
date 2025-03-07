function out=extremaldgt2ep()

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
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
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
arcarg=OCMATCONT.HE.arcarg(arc);

arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATCONT.HE.numparameter)];
diffarctime=diff(arctime);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
%dxdt=dtds(ones(OCMATCONT.DOMAINDDATA(arcindex).numeq,1),:).*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
dxdt=dtds*OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
if OCMATAE.objectivevaluecalc
    dxdt(OCMATAE.objectivevaluecoord,:)=dtds.*OCMATAE.objectivefunction(t,depvar,modelpar,arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
arcarg=OCMATCONT.HE.arcarg(arc);
arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATCONT.HE.numparameter)];
diffarctime=diff(arctime);
t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=OCMATAE.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
if OCMATAE.objectivevaluecalc
    J=[J; ...
        dtds*OCMATAE.objectivefunctionjacobian(t,depvar,modelpar,arcarg)];
    J=[J zeros(OCMATCONT.numeq,OCMATAE.playernum)];
end
Jpar=zeros(OCMATCONT.numeq,OCMATCONT.HE.numparameter);
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
        Jpar(OCMATCONT.canonicalsystemcoordinate,OCMATAE.switchtimecoord(arc))=dxdt+diffarctime(arc)*(s-arc)*Jt;
        if arc>1
            Jpar(OCMATCONT.canonicalsystemcoordinate,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        end
    else
        Jpar(OCMATCONT.canonicalsystemcoordinate,OCMATAE.switchtimecoord(arc-1))=-(dxdt+diffarctime(arc)*(s-arc-1)*Jt);
        if OCMATAE.movinghorizon
            dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATAE.movinghorizoncoord)=dxdt;
        end
    end
else
    if OCMATAE.movinghorizon
        dxdt=OCMATAE.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATAE.movinghorizoncoord)=dxdt;
    end
end
%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP

switchtimes=freepar(OCMATAE.switchtimecoord);
resconnec=[];

resinit=OCMATAE.bcinitial(depvara,OCMATAE.initialcoordinate,OCMATAE.initialstate,modelpar,OCMATCONT.HE.arcarg(1));
try
    resasym=OCMATAE.bcasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint);
catch
    resasym=OCMATAE.bcasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint,modelpar,OCMATCONT.HE.arcarg(end),depvara);
end
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    resconnec=[resconnec;
        OCMATAE.reset(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii); ...
        OCMATAE.guard(depvara,depvarb,modelpar,switchtimes,OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii)];
end
res=[resinit;resconnec;resasym];

%-------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];



%----------------------------------------------------------------
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP

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
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
    [constr labelS]=OCMATAE.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    if any(strfind(OCMATAE.pathtype,'u'))
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
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff);
arctime=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
diffarctime=diff(arctime);
for arc=1:numel(OCMATCONT.HE.arcarg)
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s(leftarcindex(arc):rightarcindex(arc))+(arctime(arc)-diffarctime(arc)*(arc-1));
end
h=OCMATAE.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATAE OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
switch OCMATAE.targettype
    case 'd'
        yend=y(OCMATCONT.equilibriumcoord,end);
        equi=OCMATAE.saddlepoint;
        fprintf(1,' Distance |yhat-yend|-d: %3.7g\n',norm(equi-yend)-OCMATAE.targetvalue);
    case 'T'
        fprintf(1,' Told-T: %g\n',freepar(OCMATCONT.HE.numparameter)-OCMATAE.targetvalue);
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
        % Jacobian extended with bordering vectors v and w
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
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

failed=[];
for ii=id
    switch ii
        case 1
            switch OCMATAE.targettype
                case 'd'
                    yend=y(OCMATCONT.equilibriumcoord,end);
                    equi=OCMATAE.saddlepoint;
                    out=norm(equi-yend)-OCMATAE.targetvalue;
                case 'T'
                    out=freepar(OCMATCONT.HE.numparameter)-OCMATAE.targetvalue;
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
global OCMATCONT OCMATAE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,OCMATAE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' freepar(OCMATCONT.HE.numparameter)];%[OCMATAE.initialtime freepar(OCMATAE.switchtimecoord).' OCMATAE.truncationtime];
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATAE.initialtime;
out.timehorizon=out.arcinterval(end);
out.modelparameter=OCMATAE.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremaldgt2ep';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATAE.pathtype;
out.solverinfo.switchtimecoord=OCMATAE.switchtimecoord;
if OCMATAE.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATAE.objectivevaluecoord;
else
    out.solverinfo.objectivevaluecoord=[];
end

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        s=[];
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
global OCMATCONT OCMATAE

modelpar=OCMATAE.parametervalue;
y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
z=[];
freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATAE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATAE=OCMATAE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATAE.basicglobalvarfilename '4extremaldgt2ep'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4extremaldgt2ep'],'sout','bvpout')
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

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;
