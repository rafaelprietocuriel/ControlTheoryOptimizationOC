function out=heteroclinic()

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
global OCMATHET OCMATCONT OCBVP
solutionindex=OCMATHET.solutionindex(arc);
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if OCMATHET.movinghorizon(solutionindex)
    arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.movinghorizoncoord(solutionindex))];
else
    if strcmp(OCMATHET.limitsettype{solutionindex},'e')
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    else
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.truncationtimecoord{solutionindex})];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
solutionindex=OCMATHET.solutionindex(arc);

if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if OCMATHET.movinghorizon(solutionindex)
    arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.movinghorizoncoord(solutionindex))];
else
    if strcmp(OCMATHET.limitsettype{solutionindex},'e')
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    else
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.truncationtimecoord{solutionindex})];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATHET.canonicalsystemjacobian(t,depvar,modelpar,arcarg);
J=dtds*J;
Jpar=zeros(OCMATCONT.DOMAINDDATA(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATHET.numarc(solutionindex)>1
    dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
    if ~OCMATHET.autonomous
        Jt=OCMATHET.canonicalsystemderivativetime(t,depvar,modelpar,arcarg);
    else
        Jt=0;
    end
    if relarc<OCMATHET.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoord{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-relarc)*Jt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        if OCMATHET.movinghorizon(solutionindex)
            dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATHET.movinghorizoncoord(solutionindex))=dxdt;
        end
    end
else
    if OCMATHET.movinghorizon(solutionindex)
        dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATHET.movinghorizoncoord(solutionindex))=dxdt;
    end
end
Jmodelpar=dtds*OCMATHET.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATHET.parametervaluecoord)=Jmodelpar(:,OCMATHET.parameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
resinit=[];
resasym=[];
resconnec=[];
resequilibrium=[];
resricatti=[];
userbc=[];
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

for ii=1:OCMATHET.hetorder
    switchtimes=freepar(OCMATHET.switchtimecoord{ii});
    arcarg=OCMATCONT.HE.arcarg(OCMATHET.arccoord{ii});
    if ii<OCMATHET.hetorder
        if length(OCMATHET.parameterindex)==1
            resinit=[resinit; ...
                depvara(OCMATHET.statecoord,OCMATHET.initcoord(ii+1))-depvara(OCMATHET.statecoord,OCMATHET.initcoord(ii)); ...
                depvara(OCMATHET.costatecoord,OCMATHET.initcoord(ii+1))-depvara(OCMATHET.costatecoord,OCMATHET.initcoord(ii))-freepar(end)*OCMATHET.initialcostatedifference];
        else
            resinit=[resinit; ...
                depvara(:,OCMATHET.initcoord(ii+1))-depvara(:,OCMATHET.initcoord(ii))];
        end
        if ~isempty(OCMATHET.fixcoordinate)
            resinit=[resinit; ...
                OCMATHET.fixvalue-depvara(OCMATHET.fixcoordinate,OCMATHET.initcoord(ii))];
        end
    end
    hatx=freepar(OCMATHET.equilibriumcoord{ii});
    resequilibrium=[resequilibrium;... 
        OCMATHET.equilibrium(hatx,modelpar,arcarg(end))];
    if OCMATHET.simple
        asymptoticmatrix=OCMATHET.asymptoticmatrix{ii};
    else
        Y=freepar(OCMATHET.Ycoord{ii});
        asymptoticmatrix=OCMATHET.Q0{ii}*[-Y';OCMATHET.Id{ii}];
        Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
        resricatti=[resricatti;ricatti(Y,Jac,OCMATHET.Q0{ii},OCMATHET.subspacedim{ii})];
    end
    resasym=[resasym;OCMATHET.bcasymptotic(depvarb(:,OCMATHET.cumsumnumarc(ii)),asymptoticmatrix,hatx)];
    for arc=1:OCMATHET.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATHET.reset(depvara(:,OCMATHET.arccoord{ii}),depvarb(:,OCMATHET.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc); ...
            OCMATHET.guard(depvara(:,OCMATHET.arccoord{ii}),depvarb(:,OCMATHET.arccoord{ii}),modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc)];
    end
    if OCMATHET.movinghorizon(ii)
        resasym=[resasym; ...
            sqrt(sum((hatx-depvarb(:,OCMATHET.cumsumnumarc(ii))).^2))-OCMATHET.distance(ii)];
    end
end
if ~isempty(OCMATHET.userbc)
    userbc=OCMATHET.userbc(depvara,depvarb,freepar,modelpar,OCMATHET);
end

res=[resinit;resasym;resconnec;resequilibrium;resricatti;userbc];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT OCBVP
Ja=[];
Jb=[];
Jpar=[];

%-------------------------------------------------------------------------
function failed=dataadaptation(tmesh,coeff,tangent,tmeshold)
failed=0;

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATHET OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counter=0;
for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATHET.solutionindex(arc);
    if solutionindex>1
        relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    if OCMATHET.movinghorizon(solutionindex)
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.movinghorizoncoord(solutionindex))];
    else
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATHET.testadmissibility(t,sol.y(idx,:),modelpar,arcarg);
    else
        eqcoord=domainddata(arcindex).eqcoord;
        [constr labelS]=OCMATHET.testadmissibility(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
    if OCMATHET.stableflag{solutionindex}
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
for ii=1:OCMATHET.hetorder
    hatx=freepar(OCMATHET.equilibriumcoord{ii});
    violationmat=OCMATCONT.OPTIONS.maxdistance-norm(y(:,rightarcindex(OCMATHET.cumsumnumarc(ii)))-hatx)<0;
    if ~OCMATHET.movinghorizon(ii)
        if violationmat
            counter=counter+1;
            cols=size(y,2);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='maxdistance';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=norm(y(:,rightarcindex(OCMATHET.cumsumnumarc(ii)))-hatx);
            infoS(counter).minval=OCMATCONT.OPTIONS.maxdistance-norm(y(:,rightarcindex(OCMATHET.cumsumnumarc(ii)))-hatx);
            b=min([b infoS(counter).minval]);
        end
    end
    if OCMATHET.testhopf
        Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
        eigval=eig(Jac);
        violationmat=max(abs(imag(eigval)))>0;
        if violationmat
            counter=counter+1;
            cols=size(y,2);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='imag_eigval';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=max(abs(imag(eigval)));
            infoS(counter).minval=max(abs(imag(eigval)));
            b=min([b -infoS(counter).minval]);
        end
    end
    if OCMATHET.testeig
        if ~OCMATHET.testhopf
            Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
            eigval=eig(Jac);
        end
        if OCMATHET.stableflag{ii}
            violationmat=abs(OCMATHET.subspacedim{ii}-length(find(eigval<0)));
        else
            violationmat=abs(OCMATHET.subspacedim{ii}-length(find(eigval>0)));
        end
        if violationmat
            counter=counter+1;
            cols=size(y,2);
            infoS(counter).arcarg=arcarg;
            infoS(counter).arcnum=arc;
            infoS(counter).rows='num_eigval';
            infoS(counter).cols=cols;
            infoS(counter).violationmat=violationmat;
            infoS(counter).constraintvalue=[];
            infoS(counter).minval=[];
            b=min([b -violationmat]);
        end
    end
end
%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);
% clear possible persistent variable
h=OCMATHET.plotcontinuation(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
set(h(OCMATHET.arccoord{2}),'LineStyle','--')
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
idx=[];
if isempty(coeff)
    return
end

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
if isempty(OCMATHET.targetparametervalue) && isempty(OCMATHET.targetfunction)
    fprintf(1,' Continuation parameter: %g\n',coeff(end));
elseif isempty(OCMATHET.targetfunction)
    fprintf(1,' Difference to targetvalue: %g\n',OCMATHET.targetparametervalue-modelpar(OCMATHET.targetparameterindex));
    fprintf(1,' Parameter: %g, %g\n',freepar(OCMATHET.parametervaluecoord(1)),freepar(OCMATHET.parametervaluecoord(2)));
else
    if ~OCMATHET.hitfunction
        out=OCMATHET.targetfunction(freepar,modelpar,OCMATHET);
    else
        depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
        depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
        out=OCMATHET.targetfunction(depvara,depvarb,modelpar,OCMATHET.arcarg);
    end
    fprintf(1,' Targetvalue: %g\n',out);
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET

failed=[];
for ii=id
    switch ii
        case 1
            if isempty(OCMATHET.targetparametervalue) && isempty(OCMATHET.targetfunction)
                out=coeff(end);
            elseif isempty(OCMATHET.targetfunction)
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                out=OCMATHET.targetparametervalue-modelpar(OCMATHET.targetparameterindex);
            else
                [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                if ~OCMATHET.hitfunction
                    out=OCMATHET.targetfunction(freepar,modelpar,OCMATHET);
                else
                    depvara=y(:,OCMATCONT.HE.TIMEDDATA.leftarcindex);
                    depvarb=y(:,OCMATCONT.HE.TIMEDDATA.rightarcindex);
                    out=OCMATHET.targetfunction(depvara,depvarb,modelpar,OCMATHET.arcarg);
                end
            end
            %out=tangent(end);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATHET OCBVP
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out=transform2nativematlab(tmesh,coeff,OCMATHET.parametervalue);
if OCBVP.multiarccalc
    out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
        OCMATCONT.HE.TIMEDDATA.rightarcindex];
else
    out.arcposition=OCMATCONT.HE.arcrowindex;
end
out.arcinterval=[];
for solutionindex=1:2
    %solutionindex=OCMATHET.solutionindex(arc);
    if OCMATHET.movinghorizon(solutionindex)
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' freepar(OCMATHET.movinghorizoncoord(solutionindex))];
    else
        arctime=[OCMATHET.initialtime freepar(OCMATHET.switchtimecoord{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon(solutionindex)=arctime(end);
    out.solverinfo.arcarg{solutionindex}=OCMATCONT.HE.arcarg(OCMATHET.arccoord{solutionindex});
    out.solverinfo.arcinterval{solutionindex}=arctime;
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATHET.initialtime;
%out.modelparameter=OCMATHET.parametervalue;
out.modelparameter=modelpar;
out.modelparameter(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.parameterindex=OCMATHET.parameterindex;
out.solverinfo.parametervaluecoord=OCMATHET.parametervaluecoord;
out.solverinfo.conttype='heteroclinic';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.inftimetransformation=OCMATHET.inftimetransformation;
out.solverinfo.pathtype=OCMATHET.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATHET.switchtimecoord;
out.solverinfo.truncationtimecoord=OCMATHET.truncationtimecoord;
out.solverinfo.equilibriumcoord=OCMATHET.equilibriumcoord;
out.solverinfo.solutionindex=OCMATHET.solutionindex;
if ~OCMATHET.simple
    out.solverinfo.Ycoord=OCMATHET.Ycoord;
    out.solverinfo.subspacedim=OCMATHET.subspacedim;
    out.solverinfo.orthspacedim=OCMATHET.orthspacedim;
    out.solverinfo.qbasis=OCMATHET.Q0;
end
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
global OCMATCONT OCMATHET
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
global OCMATCONT OCMATHET

%numarc=OCMATCONT.HE.numarc;
%domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATHET.parametervalue;
% if length(OCMATHET.parameterindex)==2
%     modelpar(OCMATHET.parameterindex)=coeff(OCMATCONT.HE.contparametercoord+(-1:0));
% else
%     modelpar(OCMATHET.parameterindex)=coeff(OCMATCONT.HE.contparametercoord-1);
% end
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATHET OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATHET=OCMATHET;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATHET.basicglobalvarfilename '4heteroclinic'],'MODELINFO')
    end
    save([OCMATHET.basicresultfilename '4heteroclinic'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function discretizationdata=domaindiscretization(arcarg)
global OCMATHET

discretizationdata=OCMATHET.domaindiscretization(arcarg);

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATHET

pathname=OCMATHET.datapath();


%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent]=adapt(tmesh,coeff,tangent)
global OCMATCONT OCMATHET

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);

for ii=1:OCMATHET.hetorder
    Y=freepar(OCMATHET.Ycoord{ii});
    OCMATCONT.adapted = 1;
    %
    [U,S,V]=svd(OCMATHET.Q0{ii}(:,1:OCMATHET.subspacedim{ii})+OCMATHET.Q0{ii}(:,OCMATHET.subspacedim{ii}+1:end)*Y);
    OCMATHET.Q0{ii}= U;
    OCMATHET.Y{ii}=zeros(OCMATHET.orthspacedim{ii},OCMATHET.subspacedim{ii});

    freepar(OCMATHET.Ycoord{ii})=OCMATHET.Y{ii};
    switch OCMATCONT.bvpmethod
        case {'bvp6c','bvp4c'}
            coeff=[y(:);freepar];
        otherwise
    end
end
flag = 1;

%-----------------------------------------------------------------
function [tmeshnew,coeffnew,tangentnew]=deval(tmesh,coeff,tangent,tmeshnew)
global OCMATCONT OCMATHET
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        tmeshnew=tmesh;
        coeffnew=coeff;
    case {'bvp4c','bvp6c'}
        [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,OCMATHET.parametervalue);
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
