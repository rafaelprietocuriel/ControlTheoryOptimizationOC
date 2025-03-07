function out=indifferencesolution4fimp()

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
global IOCMATINDIF OCMATCONT OCBVP
solutionindex=IOCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-IOCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
arctime=IOCMATINDIF.arctime{solutionindex};
arctime(IOCMATINDIF.switchtimeidx{solutionindex})=freepar(IOCMATINDIF.switchtimecoord{solutionindex});
arctime(end)=freepar(IOCMATINDIF.endtimecoord);
diffarctime=diff(arctime);
jumparg=IOCMATINDIF.jumparg(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
if solutionindex>1
    transformedtimeshift=IOCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*IOCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg,jumparg);
dxdt(IOCMATINDIF.objectivevaluecoord,:)=dtds*IOCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg,jumparg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global IOCMATINDIF OCMATCONT OCBVP
solutionindex=IOCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-IOCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
arctime=IOCMATINDIF.arctime{solutionindex};
arctime(IOCMATINDIF.switchtimeidx{solutionindex})=freepar(IOCMATINDIF.switchtimecoord{solutionindex});
arctime(end)=freepar(IOCMATINDIF.endtimecoord);
diffarctime=diff(arctime);
jumparg=IOCMATINDIF.jumparg(arc);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
if solutionindex>1
    transformedtimeshift=IOCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=IOCMATINDIF.canonicalsystemjacobian(t,depvar,modelpar,arcarg,jumparg);
J=dtds*J;
J=[J; ...
    dtds*IOCMATINDIF.objectivefunctionjacobian(t,depvar,modelpar,arcarg,jumparg)];
J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if IOCMATINDIF.varyarcintervalidx{solutionindex}(relarc) || IOCMATINDIF.varyarcintervalidx{solutionindex}(relarc+1)
    dxdt=IOCMATINDIF.canonicalsystem(t,depvar,modelpar,arcarg,jumparg);
    if ~IOCMATINDIF.autonomous
        Jt=IOCMATINDIF.canonicalsystemderivativetime(t,depvar,modelpar,arcarg,jumparg);
    else
        Jt=0;
    end
    dxdt(IOCMATINDIF.objectivevaluecoord,:)=IOCMATINDIF.objectivefunction(t,depvar,modelpar,arcarg,jumparg);
    Jt(IOCMATINDIF.objectivevaluecoord,:)=IOCMATINDIF.objectivefunctionderivativetime(t,depvar,modelpar,arcarg,jumparg);
end
if IOCMATINDIF.varyarcintervalidx{solutionindex}(relarc)
    Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,IOCMATINDIF.varyarcintervalidx{solutionindex}(relarc))=-(dxdt+diffarctime(relarc)*(s-transformedtimeshift-relarc)*Jt);
end
if IOCMATINDIF.varyarcintervalidx{solutionindex}(relarc+1)
    Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,IOCMATINDIF.varyarcintervalidx{solutionindex}(relarc+1))=dxdt+diffarctime(relarc)*(s-transformedtimeshift-relarc+1)*Jt;
end


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global IOCMATINDIF OCMATCONT
resconnec=[];
residpt=[];
resjump=[];
resinit=[];
restrans=[];
restarget=[];

O=zeros(1,IOCMATINDIF.indifferenceorder);
endtime=freepar(IOCMATINDIF.endtimecoord);
for ii=1:IOCMATINDIF.indifferenceorder
    jumparg=IOCMATINDIF.jumparg(IOCMATINDIF.jumpcoord{ii});
    arctime=IOCMATINDIF.arctime{ii};
    arctime(IOCMATINDIF.switchtimeidx{ii})=freepar(IOCMATINDIF.switchtimecoord{ii});
    arctime(end)=endtime;
    X0=freepar(IOCMATINDIF.initialdepvarcoord{ii});
    XT=freepar(IOCMATINDIF.enddepvarcoord{ii});
    if ii==1
        if ~isempty(IOCMATINDIF.freevectorcoord)
            restarget=X0(IOCMATINDIF.statecoordinate)-(IOCMATINDIF.initialstate+freepar(IOCMATINDIF.freevectorcoord)*IOCMATINDIF.freevector);
        end
    end
    arcarg=OCMATCONT.HE.arcarg(IOCMATINDIF.arccoord{ii});
    resjump=[resjump; ...
        IOCMATINDIF.bcevent(0,[X0 depvara(IOCMATINDIF.statecostatecoord,IOCMATINDIF.arccoord{ii}(1))],modelpar,jumparg(1),endtime)];
    resjump=[resjump; ...
        IOCMATINDIF.bcevent(endtime,[depvarb(IOCMATINDIF.statecostatecoord,IOCMATINDIF.arccoord{ii}(end)) XT],modelpar,jumparg(end),endtime)];
    if ii<IOCMATINDIF.indifferenceorder
        resinit=[resinit; ...
            depvara(IOCMATINDIF.statecoordinate,IOCMATINDIF.initcoord(ii+1))-depvara(IOCMATINDIF.statecoordinate,IOCMATINDIF.initcoord(ii))];
%         residpt=[residpt; ...
%             IOCMATINDIF.bcindifference(depvara,modelpar,OCMATCONT.HE.arcarg([IOCMATINDIF.arccoord{:}]),IOCMATINDIF.initcoord([ii ii+1]))];
    end
    resinit=[resinit; ...
        depvara(IOCMATINDIF.objectivevaluecoord,IOCMATINDIF.arccoord{ii}(1))];
    restrans=[restrans; ...
        IOCMATINDIF.bctransversality(endtime,[depvarb(IOCMATINDIF.statecostatecoord,IOCMATINDIF.arccoord{ii}(end)) XT],modelpar,arcarg(end),jumparg(end))];

    for arc=1:IOCMATINDIF.numarc(ii)-1
        if jumparg(arc+1)
            resconnec=[resconnec;
                IOCMATINDIF.bcevent(arctime(arc+1),[depvarb(:,IOCMATINDIF.arccoord{ii}(arc)) depvara(:,IOCMATINDIF.arccoord{ii}(arc+1))],modelpar,jumparg(arc+1),freepar(end)); ...
                IOCMATINDIF.bcinteriorevent(arctime(arc+1),[depvarb(:,IOCMATINDIF.arccoord{ii}(arc)) depvara(:,IOCMATINDIF.arccoord{ii}(arc+1))],modelpar,arcarg(arc:arc+1),jumparg(arc+1))];
        else
            resconnec=[resconnec;
                IOCMATINDIF.reset(depvara(:,IOCMATINDIF.arccoord{ii}),depvarb(:,IOCMATINDIF.arccoord{ii}),modelpar,arctime(2:end),arcarg,IOCMATINDIF.edge{ii},arc,jumparg); ...
                IOCMATINDIF.guard(depvara(:,IOCMATINDIF.arccoord{ii}),depvarb(:,IOCMATINDIF.arccoord{ii}),modelpar,arctime(2:end),arcarg,IOCMATINDIF.edge{ii},arc,jumparg)];
        end
        resconnec=[resconnec;depvarb(IOCMATINDIF.objectivevaluecoord,IOCMATINDIF.arccoord{ii}(arc))-depvara(IOCMATINDIF.objectivevaluecoord,IOCMATINDIF.arccoord{ii}(arc+1))];
    end
    O(ii)=IOCMATINDIF.objectivevalue(depvara(:,IOCMATINDIF.arccoord{ii}),depvarb(:,IOCMATINDIF.arccoord{ii}),X0,XT,modelpar,arctime,arcarg,jumparg);
end
for ii=1:IOCMATINDIF.indifferenceorder-1
         residpt=[residpt; ...
             O(ii)-O(ii+1)];
end
res=[resinit;restarget;restrans;resconnec;resjump;residpt];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT IOCMATINDIF OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=IOCMATINDIF.solutionindex(arc);
    if solutionindex>1
        relarc=arc-IOCMATINDIF.cumsumnumarc(solutionindex-1);
        transformedtimeshift=IOCMATINDIF.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    arctime=IOCMATINDIF.arctime{solutionindex};
    arctime(IOCMATINDIF.switchtimeidx{solutionindex})=freepar(IOCMATINDIF.switchtimecoord{solutionindex});
    arctime(end)=freepar(IOCMATINDIF.endtimecoord);
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=IOCMATINDIF.testadmissibility(t,sol.y(idx,:),modelpar,arcarg,[]);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=IOCMATINDIF.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg,[]);
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
global IOCMATINDIF OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
if length(sol.x)~=size(sol.y,2)
    sol
end
% clear possible persistent variable
h=IOCMATINDIF.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
%global IOCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
fprintf(1,' Continuation parameter: %g\n',coeff(end));

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT IOCMATINDIF
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
global OCMATCONT IOCMATINDIF
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
global OCMATCONT IOCMATINDIF

failed=[];
out=[];
for ii=id
    switch ii
        case 1
            if IOCMATINDIF.targettype==1
                out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
            elseif ~isempty(IOCMATINDIF.targetvalue)
                out=IOCMATINDIF.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
            end
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT IOCMATINDIF OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,IOCMATINDIF.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[];
out.timehorizon=[];
for arc=1:IOCMATINDIF.indifferenceorder
    arctime=IOCMATINDIF.arctime{arc};
    arctime(IOCMATINDIF.switchtimeidx{arc})=freepar(IOCMATINDIF.switchtimecoord{arc});
    arctime(end)=freepar(IOCMATINDIF.endtimecoord);
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=IOCMATINDIF.initialtime;
out.modelparameter=IOCMATINDIF.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution4fimp';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=IOCMATINDIF.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=IOCMATINDIF.switchtimecoord;
out.solverinfo.statecostatecoord=IOCMATINDIF.statecostatecoord;
out.solverinfo.initialdepvarcoord=IOCMATINDIF.initialdepvarcoord;
out.solverinfo.enddepvarcoord=IOCMATINDIF.enddepvarcoord;
out.solverinfo.switchtimeidx=IOCMATINDIF.switchtimeidx;
out.solverinfo.varyarcintervalidx=IOCMATINDIF.varyarcintervalidx;
out.solverinfo.numarc=IOCMATINDIF.numarc;
out.solverinfo.initialstateindex=IOCMATINDIF.initialstateindex;
out.solverinfo.jumparg=IOCMATINDIF.jumparg;
out.solverinfo.jumpcoord=IOCMATINDIF.jumpcoord;
out.solverinfo.arccoord=IOCMATINDIF.arccoord;
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
global OCMATCONT IOCMATINDIF
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        IOCMATINDIF.PD_phi=p'/norm(p);
        IOCMATINDIF.PD_psi=Q(:,end);
        s.data.phi=IOCMATINDIF.PD_phi(:);
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
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT IOCMATINDIF

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=IOCMATINDIF.parametervalue;
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

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT IOCMATINDIF OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.IOCMATINDIF=IOCMATINDIF;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([IOCMATINDIF.basicglobalvarfilename '4indifferencesolution4fimp'],'MODELINFO')
    end
    save([IOCMATINDIF.basicresultfilename '4indifferencesolution4fimp'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global IOCMATINDIF

discretizationdata=IOCMATINDIF.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global IOCMATINDIF

pathname=IOCMATINDIF.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT IOCMATINDIF
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,IOCMATINDIF.parametervalue);
        nN=n*N;
        % under the giving assumptions the solution is continuous even at
        % the border of the arcs. Therefore, the correspondent warning of
        % deval is turned 'off'
        mbcidx=find(~diff(tmesh));
        Lidx = [1, mbcidx+1];
        Ridx = [mbcidx, length(tmesh)];
        mbcidx=find(~diff(tmeshnew));
        Lidxnew = [1, mbcidx+1];
        Ridxnew = [mbcidx, length(tmeshnew)];
        warning('off','MATLAB:deval:NonuniqueSolution');
        ynew=zeros(size(sol.y,1),length(tmeshnew));
        solpart.solver=sol.solver;
        for ii=1:length(Lidx)
            solpart.y=sol.y(:,Lidx(ii):Ridx(ii));
            solpart.yp=sol.yp(:,Lidx(ii):Ridx(ii));
            solpart.x=sol.x(Lidx(ii):Ridx(ii));
            ynew(:,Lidxnew(ii):Ridxnew(ii))=devalbvpoc(solpart,tmeshnew(Lidxnew(ii):Ridxnew(ii)));
        end
        warning('on','MATLAB:deval:NonuniqueSolution');
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,n,N,nN);

        coeffnew=[ynew(:);freepar(:)];
    case 'bvp5c'
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,IOCMATINDIF.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end