function out=extremalt4fimp()

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

function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function dxdt=ode(s,depvar,arc,freepar,modelpar)
global IOCMATFTE OCMATCONT OCBVP
arctime=IOCMATFTE.arctime;
arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);

diffarctime=diff(arctime);
jumparg=IOCMATFTE.jumparg(arc);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
dxdt=dtds*IOCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg,jumparg);
if IOCMATFTE.objectivevaluecalc
    dxdt(IOCMATFTE.objectivevaluecoord,:)=dtds*IOCMATFTE.objectivefunction(t,depvar,modelpar,arcarg,jumparg);
end
% if  any(abs(imag(dxdt))>0)
%     dxdt
% end


%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global IOCMATFTE OCMATCONT OCBVP
arctime=IOCMATFTE.arctime;
arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);

diffarctime=diff(arctime);
arcindex=OCMATCONT.HE.arcindex(arc);
jumparg=IOCMATFTE.jumparg(arc);
arcarg=OCMATCONT.HE.arcarg(arc);

t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1));
dtds=diffarctime(arc);
J=IOCMATFTE.canonicalsystemjacobian(t,depvar,modelpar,arcarg,jumparg);
J=dtds*J;
if IOCMATFTE.objectivevaluecalc
    J=[J; ...
        dtds*IOCMATFTE.objectivefunctionjacobian(t,depvar,modelpar,arcarg,jumparg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if IOCMATFTE.varyarcintervalidx(arc) || IOCMATFTE.varyarcintervalidx(arc+1)
    dxdt=IOCMATFTE.canonicalsystem(t,depvar,modelpar,arcarg,jumparg);
    if ~IOCMATFTE.autonomous
        Jt=IOCMATFTE.canonicalsystemderivativetime(t,depvar,modelpar,arcarg,jumparg);
    else
        Jt=0;
    end
    if IOCMATFTE.objectivevaluecalc
        dxdt(IOCMATFTE.objectivevaluecoord,:)=IOCMATFTE.objectivefunction(t,depvar,modelpar,arcarg,jumparg);
        Jt(IOCMATFTE.objectivevaluecoord,:)=IOCMATFTE.objectivefunctionderivativetime(t,depvar,modelpar,arcarg,jumparg);
    end
end
if IOCMATFTE.varyarcintervalidx(arc)
    Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,IOCMATFTE.varyarcintervalidx(arc))=-(dxdt+diffarctime(arc)*(s-arc)*Jt);
end
if IOCMATFTE.varyarcintervalidx(arc+1)
    Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,IOCMATFTE.varyarcintervalidx(arc+1))=dxdt+diffarctime(arc)*(s-arc+1)*Jt;
end

% if any(abs(imag(J(:)))>0) || any(abs(imag(Jpar(:)))>0)
%     J,Jpar
% end

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global IOCMATFTE OCMATCONT OCBVP
arctime=IOCMATFTE.arctime;
arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);

X0=freepar(IOCMATFTE.initialdepvarcoord);
XT=freepar(IOCMATFTE.enddepvarcoord);

resconnec=[];

if IOCMATFTE.objectivevaluecalc
    OVal=IOCMATFTE.salvagevalue(arctime(end),[depvarb(IOCMATFTE.initialdepvarcoord,end) XT],modelpar,OCMATCONT.HE.arcarg(end),IOCMATFTE.jumparg(end));
    if IOCMATFTE.jumparg(1)
        OVal=Oval+IOCMATFTE.impulseobjectivefunction(arctime(1),[X0 depvara(IOCMATFTE.initialdepvarcoord,1)],modelpar,OCMATCONT.HE.arcarg,IOCMATFTE.jumparg(1));
    end
end
resjump=IOCMATFTE.bcevent(0,[X0 depvara(IOCMATFTE.initialdepvarcoord,1)],modelpar,IOCMATFTE.jumparg(1),freepar(end));
resjump=[resjump; ...
    IOCMATFTE.bcevent(arctime(end),[depvarb(IOCMATFTE.initialdepvarcoord,end) XT],modelpar,IOCMATFTE.jumparg(end),freepar(end))];

resinit=IOCMATFTE.bcinitial(X0,IOCMATFTE.initialcoordinate,IOCMATFTE.initialstate,modelpar,OCMATCONT.HE.arcarg(1));
restrans=IOCMATFTE.bctransversality(arctime(end),[depvarb(IOCMATFTE.initialdepvarcoord,end) XT],modelpar,OCMATCONT.HE.arcarg(end),IOCMATFTE.jumparg(end));
for ii=1:numel(OCMATCONT.HE.arcarg)-1
    if IOCMATFTE.jumparg(ii+1)
        resconnec=[resconnec;
            IOCMATFTE.bcevent(arctime(ii+1),[depvarb(:,ii) depvara(:,ii+1)],modelpar,IOCMATFTE.jumparg(ii+1),freepar(end)); ...
            IOCMATFTE.bcinteriorevent(arctime(ii+1),[depvarb(:,ii) depvara(:,ii+1)],modelpar,OCMATCONT.HE.arcarg(ii:ii+1),IOCMATFTE.jumparg(ii+1))];
        if IOCMATFTE.objectivevaluecalc
            OVal=OVal+IOCMATFTE.impulseobjectivefunction(arctime(ii+1),[depvarb(IOCMATFTE.initialdepvarcoord,ii) depvara(IOCMATFTE.initialdepvarcoord,ii+1)],modelpar,OCMATCONT.HE.arcarg,IOCMATFTE.jumparg(ii+1));
        end
    else
        resconnec=[resconnec;
            IOCMATFTE.reset(depvara,depvarb,modelpar,arctime(2:end),OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii,IOCMATFTE.jumparg); ...
            IOCMATFTE.guard(depvara,depvarb,modelpar,arctime(2:end),OCMATCONT.HE.arcarg,OCMATCONT.HE.edge,ii,IOCMATFTE.jumparg)];
    end
    if IOCMATFTE.objectivevaluecalc
        resconnec=[resconnec;depvarb(IOCMATFTE.objectivevaluecoord,ii)-depvara(IOCMATFTE.objectivevaluecoord,ii+1)];
    end
end
if IOCMATFTE.objectivevaluecalc
    OVal=OVal+IOCMATFTE.impulseobjectivefunction(arctime(end),[depvarb(IOCMATFTE.initialdepvarcoord,end) XT],modelpar,OCMATCONT.HE.arcarg,IOCMATFTE.jumparg(end));
    resinit=[resinit;depvara(IOCMATFTE.objectivevaluecoord,1)-OVal];
end

res=[resinit;resconnec;restrans;resjump];
% if  any(abs(imag(res))>0)
%     res
% end

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global IOCMATFTE OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT IOCMATFTE OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
arctime=IOCMATFTE.arctime;
arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);

% if IOCMATFTE.varyendtime
%     arctime=[IOCMATFTE.initialtime freepar(IOCMATFTE.switchtimecoord).'];
% else
%     arctime=[IOCMATFTE.initialtime freepar(IOCMATFTE.switchtimecoord).' IOCMATFTE.arctime(end)];
% end
for arc=1:OCMATCONT.HE.numarc
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(arc)*s+(arctime(arc)-diffarctime(arc)*(arc-1)); 
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=IOCMATFTE.testadmissibility(t,sol.y(idx,:),modelpar,arcarg,[]);
    else
        %eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=IOCMATFTE.testadmissibility(t,sol.y(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg,[]);
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
function h=plotcontinuation(tmesh,coeff,tangent)
global IOCMATFTE OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
sol=evalatmesh(tmesh,y,z);

arctime=IOCMATFTE.arctime;
arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);
diffarctime=diff(arctime);

t=[];
for ii=1:length(arctime)-1
    t=[t diffarctime(ii)*sol.x(OCMATCONT.HE.TIMEDDATA.leftarcindex(ii):OCMATCONT.HE.TIMEDDATA.rightarcindex(ii))+(arctime(ii)-diffarctime(ii)*(ii-1))];
end
% clear possible persistent variable
figure(1)
h=IOCMATFTE.plotcontinuation(t,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT IOCMATFTE
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar,ya,yb]=drearr(tmesh,coeff);
fprintf(1,' Continuation parameter: %g\n',coeff(OCMATCONT.HE.contparametercoord));
if IOCMATFTE.findoptimaljumptime
    arctime=IOCMATFTE.arctime;
    arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);
    X0=freepar(IOCMATFTE.initialdepvarcoord);
    XT=freepar(IOCMATFTE.enddepvarcoord);
    arcarg=OCMATCONT.HE.arcarg([1 1:end end]);
    if IOCMATFTE.objectivevaluecalc
        ya=[ya(1:end-1,:) XT]; % right side limits
        yb=[X0 yb(1:end-1,:)]; % left side limits
    else
        ya=[ya XT]; % right side limits
        yb=[X0 yb]; % left side limits
    end
    val=IOCMATFTE.bcinteriorevent(arctime(IOCMATFTE.findoptimaljumptime),[yb(:,IOCMATFTE.findoptimaljumptime) ya(:,IOCMATFTE.findoptimaljumptime)],modelpar,arcarg([IOCMATFTE.findoptimaljumptime IOCMATFTE.findoptimaljumptime+1]),abs(IOCMATFTE.jumparg(IOCMATFTE.findoptimaljumptime)));
    fprintf(1,' Optimal Time condition: %g\n',val);
end
%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT IOCMATFTE
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
global OCMATCONT IOCMATFTE
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
global OCMATCONT IOCMATFTE 
[t,y,z,freepar,modelpar,ya,yb]=drearr(tmesh,coeff);

failed=[];
for ii=id
    switch ii
        case 1

            out=freepar(OCMATCONT.HE.numparameter)-IOCMATFTE.targetvalue;
        case 2
            arctime=IOCMATFTE.arctime;
            arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);
            X0=freepar(IOCMATFTE.initialdepvarcoord);
            XT=freepar(IOCMATFTE.enddepvarcoord);
            arcarg=OCMATCONT.HE.arcarg([1 1:end end]);
            if IOCMATFTE.objectivevaluecalc
                ya=[ya(1:end-1,:) XT]; % right side limits
                yb=[X0 yb(1:end-1,:)]; % left side limits
            else
                ya=[ya XT]; % right side limits
                yb=[X0 yb]; % left side limits
            end
            out(2)=IOCMATFTE.bcinteriorevent(arctime(IOCMATFTE.findoptimaljumptime),[yb(:,IOCMATFTE.findoptimaljumptime) ya(:,IOCMATFTE.findoptimaljumptime)],modelpar,arcarg([IOCMATFTE.findoptimaljumptime IOCMATFTE.findoptimaljumptime+1]),abs(IOCMATFTE.jumparg(IOCMATFTE.findoptimaljumptime)));
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT IOCMATFTE OCBVP
dataadaptation(tmesh);
[t,y,z,freepar,modelpar,ya,yb]=drearr(tmesh,coeff);

out=transform2nativematlab(tmesh,coeff,IOCMATFTE.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
% adaptation for hybrid solution structure
out.x=[out.x(1) out.x out.x(end)];
if IOCMATFTE.objectivevaluecalc
    out.y=[[freepar(IOCMATFTE.initialdepvarcoord);0] out.y [freepar(IOCMATFTE.enddepvarcoord);out.y(end,end)]];
else
    out.y=[freepar(IOCMATFTE.initialdepvarcoord) out.y freepar(IOCMATFTE.enddepvarcoord)];
end
out.arcposition=1+out.arcposition;
arctime=IOCMATFTE.arctime;
arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);

% if IOCMATFTE.varyendtime
%     arctime=[IOCMATFTE.initialtime freepar(IOCMATFTE.switchtimecoord).'];
% else
%     arctime=[IOCMATFTE.initialtime freepar(IOCMATFTE.switchtimecoord).' IOCMATFTE.endtime];
% end
out.arcinterval=arctime;
%out.arcinterval=[IOCMATFTE.initialtime freepar(IOCMATFTE.switchtimecoord).' IOCMATFTE.truncationtime];
out.arcarg=OCMATCONT.HE.arcarg;
out.jumparg=IOCMATFTE.jumparg;
out.x0=IOCMATFTE.initialtime;
out.timehorizon=arctime(end);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='extremalt4fimp';
out.solverinfo.objectivevaluecalc=IOCMATFTE.objectivevaluecalc;
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=IOCMATFTE.switchtimecoord;
out.solverinfo.switchtimeidx=IOCMATFTE.switchtimeidx;
out.solverinfo.initialdepvarcoord=IOCMATFTE.initialdepvarcoord;
out.solverinfo.enddepvarcoord=IOCMATFTE.enddepvarcoord;
out.solverinfo.objectivevaluecoord=IOCMATFTE.objectivevaluecoord;
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        out.solverinfo.yp=out.yp;
        out.solverinfo.ypmid=out.ypmid;
        out=rmfield(out,{'yp','ypmid'});
    case 'bvp4c'
        out.solverinfo.yp=out.yp;
        out=rmfield(out,'yp');
end
%if IOCMATFTE.findoptimaljumptime
    arctime=IOCMATFTE.arctime;
    arctime(IOCMATFTE.switchtimeidx)=freepar(IOCMATFTE.switchtimecoord);
    X0=freepar(IOCMATFTE.initialdepvarcoord);
    XT=freepar(IOCMATFTE.enddepvarcoord);
    arcarg=OCMATCONT.HE.arcarg([1 1:end end]);
    if IOCMATFTE.objectivevaluecalc
        ya=[ya(1:end-1,:) XT]; % right side limits
        yb=[X0 yb(1:end-1,:)]; % left side limits
    else
        ya=[ya XT]; % right side limits
        yb=[X0 yb]; % left side limits
    end
    if IOCMATFTE.findoptimaljumptime
        val=IOCMATFTE.bcinteriorevent(arctime(IOCMATFTE.findoptimaljumptime),[yb(:,IOCMATFTE.findoptimaljumptime) ya(:,IOCMATFTE.findoptimaljumptime)],modelpar,arcarg([IOCMATFTE.findoptimaljumptime IOCMATFTE.findoptimaljumptime+1]),abs(IOCMATFTE.jumparg(IOCMATFTE.findoptimaljumptime)));
    else
        val=IOCMATFTE.bcinteriorevent(arctime(IOCMATFTE.varyparameterindex),[yb(:,IOCMATFTE.varyparameterindex) ya(:,IOCMATFTE.varyparameterindex)],modelpar,arcarg([IOCMATFTE.varyparameterindex IOCMATFTE.varyparameterindex+1]),abs(IOCMATFTE.jumparg(IOCMATFTE.varyparameterindex)));
    end
    out.solverinfo.bcinteriorevent=val;
%end
%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT IOCMATFTE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        IOCMATFTE.PD_phi=p'/norm(p);
        IOCMATFTE.PD_psi=Q(:,end);
        s.data.phi=IOCMATFTE.PD_phi(:);
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
function [tmesh,y,z,freepar,modelpar,ya,yb]=drearr(tmesh,coeff)
global OCMATCONT IOCMATFTE OCBVP

modelpar=IOCMATFTE.parametervalue;
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        numarc=OCMATCONT.HE.numarc;
        domainddata=OCMATCONT.DOMAINDDATA;
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
        ya=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
        yb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        % handle case of pure state constraints
    otherwise
end


%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT IOCMATFTE OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.IOCMATFTE=IOCMATFTE;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([IOCMATFTE.basicglobalvarfilename '4extremalt4fimp'],'MODELINFO')
    end
    save([IOCMATFTE.basicresultfilename '4extremalt4fimp'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global IOCMATFTE

pathname=IOCMATFTE.datapath();

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global IOCMATFTE

discretizationdata=IOCMATFTE.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT IOCMATFTE
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,IOCMATFTE.parametervalue);
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
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,IOCMATFTE.parametervalue);
        neqn=n+OCMATCONT.HE.numparameter;
        neqnN=neqn*N;
        sol.idata.yp=sol.yp;
        sol.idata.ymid=sol.ymid;
        ynew=[devalbvpoc(sol,tmeshnew); ...
            freepar(:,ones(1,numel(tmeshnew)))];
        tangentnew=devaltangent(tmeshnew,tmesh,tangent,neqn,N,neqnN);
        coeffnew=ynew(:);
end