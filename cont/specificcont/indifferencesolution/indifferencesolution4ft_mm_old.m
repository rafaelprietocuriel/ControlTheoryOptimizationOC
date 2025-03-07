function out=indifferencesolution4ft_mm()

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
global OCMATINDIF OCMATCONT OCBVP
if OCMATINDIF.targettype==3 || OCMATINDIF.freeparameter
    for ii=1:OCMATINDIF.numberofmodels
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=freepar(OCMATINDIF.parametercoord);
    end
end
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if OCMATINDIF.targettype==2
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.endtimecoord)];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATINDIF.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
if OCMATINDIF.objectivevaluecalc
    dxdt(OCMATINDIF.objectivevaluecoord,:)=dtds*OCMATINDIF.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
end

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
solutionindex=OCMATINDIF.solutionindex(arc);
if solutionindex>1
    relarc=arc-OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    relarc=arc;
end
if OCMATINDIF.targettype==2
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.endtimecoord)];
else
    arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
end
diffarctime=diff(arctime);
arcarg=OCMATCONT.HE.arcarg(arc);
arcindex=OCMATCONT.HE.arcindex(arc);
if solutionindex>1
    transformedtimeshift=OCMATINDIF.cumsumnumarc(solutionindex-1);
else
    transformedtimeshift=0;
end
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
dtds=diffarctime(relarc(ones(1,numel(s))));
J=OCMATINDIF.canonicalsystemjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
J=dtds*J;
if OCMATINDIF.objectivevaluecalc
    J=[J; ...
        dtds*OCMATINDIF.objectivefunctionjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg)];
    J=[J zeros(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).numeq,1)];
end
Jpar=zeros(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).numeq,OCMATCONT.HE.numparameter);
if OCMATINDIF.numarc(solutionindex)>1
    dxdt=OCMATINDIF.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    if ~OCMATINDIF.autonomous{solutionindex}
        Jt=OCMATINDIF.canonicalsystemderivativetime{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    else
        Jt=0;
    end
    if OCMATINDIF.objectivevaluecalc
        dxdt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunction{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jt(OCMATINDIF.objectivevaluecoord,:)=OCMATINDIF.objectivefunctionderivativetime{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    end
    if relarc<OCMATINDIF.numarc(solutionindex)
        Jpar(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-relarc)*Jt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc-1)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA{solutionindex}(arcindex).eqcoord,OCMATINDIF.switchtimecoord{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc-1)*Jt);
        if OCMATINDIF.targettype==2
            dxdt=OCMATINDIF.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
            Jpar(:,OCMATINDIF.endtimecoord)=dxdt;
        end
    end
else
    if OCMATINDIF.targettype==2
        dxdt=OCMATINDIF.canonicalsystem{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jpar(:,OCMATINDIF.endtimecoord)=dxdt;
    end
end
if OCMATINDIF.targettype==3 || OCMATINDIF.freeparameter
    Jmodelpar=dtds*OCMATINDIF.canonicalsystemparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
    Jpar(OCMATINDIF.statecostatecoord{solutionindex},OCMATINDIF.parametercoord)=Jmodelpar(:,OCMATINDIF.parameterindex{solutionindex});
    if OCMATINDIF.objectivevaluecalc
        Jobjmodelpar=dtds*OCMATINDIF.objectivefunctionparameterjacobian{solutionindex}(t,depvar,modelpar{solutionindex},arcarg);
        Jpar(OCMATINDIF.objectivevaluecoord,OCMATINDIF.parametercoord)=Jobjmodelpar(:,OCMATINDIF.parameterindex{solutionindex});
    end
end


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT
resconnec=[];
residpt=[];
resinit=[];
restrans=[];
restarget=[];
resuser=[];
if OCMATINDIF.targettype==3 || OCMATINDIF.freeparameter
    for ii=1:OCMATINDIF.numberofmodels
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=freepar(OCMATINDIF.parametercoord);
    end
end

O=zeros(1,OCMATINDIF.indifferenceorder);
if OCMATINDIF.targettype==2
    endtime=freepar(OCMATINDIF.endtimecoord);
else
    endtime=OCMATINDIF.endtime;
end
for ii=1:OCMATINDIF.indifferenceorder
    switchtimes=freepar(OCMATINDIF.switchtimecoord{ii});
    if ii==1
        initialstate=OCMATINDIF.startvalue;
        if ~isempty(OCMATINDIF.freevectorcoord)
            for jj=1:length(OCMATINDIF.freevectorcoord)
                initialstate=initialstate+freepar(OCMATINDIF.freevectorcoord(jj))*OCMATINDIF.freevector(:,jj);
            end
        end
        restarget=depvara(OCMATINDIF.statecoordinate{ii},1)-initialstate;
    end
    arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii});
    if ii<OCMATINDIF.indifferenceorder
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate{ii},OCMATINDIF.initcoord(ii+1))-depvara(OCMATINDIF.statecoordinate{ii},OCMATINDIF.initcoord(ii))];
    end
    if OCMATINDIF.objectivevaluecalc
        OVal=OCMATINDIF.salvagevalue{ii}(endtime,depvarb(OCMATINDIF.statecostatecoord{ii},OCMATINDIF.arccoord{ii}(end)),modelpar{ii},arcarg(end));
        resinit=[resinit; ...
            depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(1))-OVal];
    else
        OVal=OCMATINDIF.salvagevalue{ii}(endtime,depvarb(OCMATINDIF.statecostatecoord,OCMATINDIF.arccoord{ii}(end)),modelpar{ii},arcarg(end));
        O(ii)=OVal+OCMATINDIF.objectivevalue{ii}([0 endtime],[depvara(:,OCMATINDIF.initcoord(ii)) depvarb(:,OCMATINDIF.endcoord(ii))],modelpar{ii},OCMATCONT.HE.arcarg([OCMATINDIF.arccoord{ii}]));
    end
    restrans=[restrans; ...
        OCMATINDIF.bctransversality{ii}(endtime,depvarb(OCMATINDIF.statecostatecoord{ii},OCMATINDIF.arccoord{ii}(end)),modelpar{ii},arcarg(end))];
    
    for arc=1:OCMATINDIF.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATINDIF.reset{ii}(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar{ii},switchtimes,arcarg,OCMATINDIF.edge{ii},arc); ...
            OCMATINDIF.guard{ii}(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar{ii},switchtimes,arcarg,OCMATINDIF.edge{ii},arc)];
        if OCMATINDIF.objectivevaluecalc
            resconnec=[resconnec; ...
                depvarb(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(arc))-depvara(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(arc+1))];
        end
    end
    if OCMATINDIF.objectivevaluecalc
        O(ii)=depvarb(OCMATINDIF.objectivevaluecoord,OCMATINDIF.arccoord{ii}(end));
    end
    if ~isempty(OCMATINDIF.userbc)
        resuser=[resuser OCMATINDIF.userbc{ii}(depvara(:,OCMATINDIF.arccoord{ii}),depvarb(:,OCMATINDIF.arccoord{ii}),modelpar{ii},arcarg,OCMATINDIF.userbcvalue);];
    end
end
for ii=1:OCMATINDIF.indifferenceorder-1
    residpt=[residpt; ...
        O(ii)-O(ii+1)];
end

res=[resinit;restarget;restrans;resconnec;residpt;resuser];
%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
Ja=[];
Jb=[];
Jpar=[];

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
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
    if OCMATINDIF.targettype==2
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' freepar(OCMATINDIF.endtimecoord)];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{solutionindex}).' OCMATINDIF.endtime];
    end
    diffarctime=diff(arctime);
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    s=sol.x(leftarcindex(arc):rightarcindex(arc));
    t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));
    if ~OCBVP.multiarccalc
        idx=OCMATCONT.HE.arcrowindex(1,arc):OCMATCONT.HE.arcrowindex(2,arc);
        [constr labelS]=OCMATINDIF.testadmissibility{solutionindex}(t,sol.y(idx,:),modelpar{solutionindex},arcarg);
    else
        eqcoord=domainddata{solutionindex}(arcindex).eqcoord;
        [constr labelS]=OCMATINDIF.testadmissibility{solutionindex}(t,sol.y(eqcoord,leftarcindex(arc):rightarcindex(arc)),modelpar{solutionindex},arcarg);
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
if length(sol.x)~=size(sol.y,2)
    sol
end
% clear possible persistent variable
h=OCMATINDIF.plotcontinuation{1}(sol.x,sol.y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
if ~isempty(OCMATINDIF.targetvalue)
    if OCMATINDIF.targettype==1
        diffval=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
    elseif OCMATINDIF.targettype==3
        if ~isempty(OCMATINDIF.targetvalue)
            diffval=OCMATINDIF.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
        else
            diffval=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
        end
    elseif OCMATINDIF.targettype==2
        diffval=OCMATINDIF.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
    end
    fprintf(1,' Difference to target value: %g\n',diffval);
elseif ~isempty(OCMATINDIF.hitvalue)
    for order=1
        diffval=OCMATINDIF.hitvaluefunc{1}(t,y,modelpar{1},OCMATCONT.HE.arcarg,OCMATINDIF.hitvalue);
    end
    fprintf(1,' Difference to target value: %g\n',diffval);
else
    fprintf(1,' Continuation parameter: %g\n',coeff(end));
end
%fprintf(1,' Continuation parameter: %g\n',coeff(end));

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
out=[];
for ii=id
    switch ii
        case 1
            if OCMATINDIF.targettype==1
                [t,y]=drearr(tmesh,coeff,tangent);
                out=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
            elseif OCMATINDIF.targettype==3
                if ~isempty(OCMATINDIF.targetvalue)
                    out=OCMATINDIF.targetvalue-coeff(OCMATCONT.HE.contparametercoord);
                elseif ~isempty(OCMATINDIF.hitvalue)
                    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                    for order=1
                            out=OCMATINDIF.hitvaluefunc{1}(t,y,modelpar{1},OCMATCONT.HE.arcarg,OCMATINDIF.hitvalue);
                    end
                else
                    [t,y]=drearr(tmesh,coeff,tangent);
                    out=OCMATINDIF.targetvalue-y(OCMATINDIF.targetcoordinate,1);
                end
            elseif ~isempty(OCMATINDIF.hitvalue)
                [t,y]=drearr(tmesh,coeff,tangent);
                    for order=1
                            out=OCMATINDIF.hitvaluefunc{1}(t,y,modelpar{1},OCMATCONT.HE.arcarg,OCMATINDIF.hitvalue);
                    end
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
out.arcinterval=[];
out.timehorizon=[];
for order=1:OCMATINDIF.indifferenceorder
    if OCMATINDIF.targettype==2
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{order}).' freepar(OCMATINDIF.endtimecoord)];
    else
        arctime=[OCMATINDIF.initialtime freepar(OCMATINDIF.switchtimecoord{order}).' OCMATINDIF.endtime];
    end
    out.arcinterval=[out.arcinterval arctime];
    out.timehorizon=[out.timehorizon arctime(end)];
end
if OCMATINDIF.targettype==3 || OCMATINDIF.freeparameter
    for ii=1:OCMATINDIF.numberofmodels
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=freepar(OCMATINDIF.parametercoord);
    end
end

out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATINDIF.initialtime;
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='indifferencesolution4ft_mm';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.initcoord=OCMATINDIF.initcoord;
out.solverinfo.solutionindex=OCMATINDIF.solutionindex;
out.solverinfo.initialstateindex=OCMATINDIF.initialstateindex;
out.solverinfo.indifferenceorder=OCMATINDIF.indifferenceorder;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.switchtimecoord=OCMATINDIF.switchtimecoord;
out.solverinfo.statecostatecoord=OCMATINDIF.statecostatecoord;
if OCMATINDIF.objectivevaluecalc
    out.solverinfo.objectivevaluecoord=OCMATINDIF.objectivevaluecoord;
end
out.solverinfo.parameterindex=OCMATINDIF.parameterindex;
if OCMATINDIF.targettype==3
    out.solverinfo.parametercoord=OCMATINDIF.parametercoord;
end
out.solverinfo.meshvalcoord=OCMATCONT.HE.DDATA.meshvalcoord;
out.solverinfo.parametercoord=OCMATCONT.HE.parametercoord;

out.solverinfo.numarc=OCMATINDIF.numarc;
out.solverinfo.initialstateindex=OCMATINDIF.initialstateindex;
out.solverinfo.arccoord=OCMATINDIF.arccoord;
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
global OCMATCONT OCMATINDIF
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATINDIF.PD_phi=p'/norm(p);
        OCMATINDIF.PD_psi=Q(:,end);
        s.data.phi=OCMATINDIF.PD_phi(:);
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
global OCMATCONT OCMATINDIF

numarc=OCMATCONT.HE.numarc;
domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATINDIF.parametervalue;
switch OCMATCONT.bvpmethod
    case {'bvp6c','bvp5c','bvp4c'}
        y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
        z=[];
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
    otherwise
end
if OCMATINDIF.targettype==3 || OCMATINDIF.freeparameter
    for ii=1:OCMATINDIF.numberofmodels
        modelpar{ii}(OCMATINDIF.parameterindex{ii})=freepar(OCMATINDIF.parametercoord);
    end
end

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATINDIF OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATINDIF=OCMATINDIF;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATINDIF.basicglobalvarfilename '4indifferencesolution4ft_mm'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4indifferencesolution4ft_mm'],'sout','bvpout')
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