function out=heteroclinicep2lp()

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
global OCMATCONT
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff);
[J,M]=calc_RHSJac4LC(t,y,z,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icfunjac);
OCMATCONT.monodromy=M;

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
limitcycleindex=OCMATHET.limitcycleindex(arc);
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if limitcycleindex
    arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(solutionindex))];
else
    if OCMATHET.freeendtime(solutionindex)
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
    else
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
    end
end
diffarctime=diff(arctime);
arcarg=OCMATHET.arcarg{solutionindex}(relarc);
t=diffarctime(relarc)*(s-transformedtimeshift)+(arctime(relarc)-diffarctime(relarc)*(relarc-1));

dtds=diffarctime(relarc);
dxdt=dtds*OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J,Jpar]=odejac(s,depvar,arc,freepar,modelpar)
global OCMATHET OCMATCONT
solutionindex=OCMATHET.solutionindex(arc);
limitcycleindex=OCMATHET.limitcycleindex(arc);
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

if solutionindex>1
    relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
    transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
else
    relarc=arc;
    transformedtimeshift=0;
end
if limitcycleindex
    arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(solutionindex))];
else
    if OCMATHET.freeendtime(solutionindex)
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
    else
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
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
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoordinate{solutionindex}(relarc))=dxdt+diffarctime(relarc)*(s-relarc)*Jt;
        if relarc>1
            Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        end
    else
        Jpar(OCMATCONT.DOMAINDDATA(arcindex).eqcoord,OCMATHET.switchtimecoordinate{solutionindex}(relarc-1))=-(dxdt+diffarctime(relarc)*(s-relarc)*Jt);
        if limitcycleindex
            dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATHET.periodcoordinate(solutionindex))=dxdt;
        else
            if OCMATHET.freeendtime(solutionindex)
                dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
                Jpar(:,OCMATHET.truncationtimecoordinate(solutionindex))=dxdt;
            end
        end
    end
else
    if limitcycleindex
        dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
        Jpar(:,OCMATHET.periodcoordinate(solutionindex))=dxdt;
    else
        if OCMATHET.freeendtime(solutionindex)
            dxdt=OCMATHET.canonicalsystem(t,depvar,modelpar,arcarg);
            Jpar(:,OCMATHET.truncationtimecoordinate(solutionindex))=dxdt;
        end
    end
end
Jmodelpar=dtds*OCMATHET.canonicalsystemparameterjacobian(t,depvar,modelpar,arcarg);
Jpar(:,OCMATHET.parametervaluecoord)=Jmodelpar(:,OCMATHET.parameterindex);


%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATHET OCMATCONT
resasym=[];
resinit=[];
resconnec=[];
resequilibrium=[];
resricatti=[];
modelpar(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);

ctr=0;
initialdepvar=zeros(OCMATHET.statecostatecoordinate(end),2);
depvaralc=depvara(:,OCMATHET.limitcycleindex==1);
depvarblc=depvarb(:,OCMATHET.limitcycleindex==1);
resper=OCMATHET.bclimitcycle(depvaralc,depvarblc);
resper=[resper; ...
    sum(OCMATHET.velocityvector.*(depvaralc(OCMATHET.velocitycoordinate,1)-OCMATHET.initialpoint))];

for ii=1:OCMATHET.hetorder
    solutionindex=find(OCMATHET.solutionindex==ii);
    limitcycleindex=OCMATHET.limitcycleindex(solutionindex(1));
    actdepvara=depvara(:,solutionindex);
    actdepvarb=depvarb(:,solutionindex);
    arcarg=OCMATHET.arcarg{ii};
    switchtimes=freepar(OCMATHET.switchtimecoordinate{ii});

    if ~limitcycleindex
        ctr=ctr+1;
        initialdepvar(:,ctr)=actdepvara(:,1);
        switch OCMATHET.limitsettype{ctr}
            case 'e'
                hatx=freepar(OCMATHET.equilibriumcoordinate{ctr});
                resequilibrium=[resequilibrium;OCMATHET.canonicalsystem(0,hatx,modelpar,arcarg(end))];
                Jac=OCMATHET.canonicalsystemjacobian(0,hatx,modelpar,arcarg(end));
            case 'l'
                Jac=OCMATCONT.monodromy;
                hatx=depvaralc(:,1);
        end
        if OCMATHET.simple
            asymptoticmatrix=OCMATHET.asymptoticmatrix{ctr};
        else
            Y=freepar(OCMATHET.Ycoordinate{ctr});
            asymptoticmatrix=OCMATHET.Q0{ctr}*[-Y';OCMATHET.Id{ctr}];
            resricatti=[resricatti;ricatti(Y,Jac,OCMATHET.Q0{ctr},OCMATHET.subspacedim{ctr})];
        end
        resasym=[resasym;OCMATHET.bcasymptotic(actdepvarb(:,end),asymptoticmatrix,hatx)];
        if OCMATHET.fixdistance(ctr)
            resasym=[resasym; ...
                sqrt(sum((hatx-actdepvarb(:,end)).^2))-OCMATHET.distance(ctr)];
        end
    end
    for arc=1:OCMATHET.numarc(ii)-1
        resconnec=[resconnec; ...
            OCMATHET.reset(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc); ...
            OCMATHET.guard(actdepvara,actdepvarb,modelpar,switchtimes,arcarg,OCMATHET.edge{ii},arc)];
    end
    if ii==1
        resinit=[resinit; ...
            actdepvara(OCMATHET.fixcoordinate,1)-OCMATHET.fixvalue];
    end
end
if OCMATHET.findinitconnection
    resinit=[resinit; ...
        initialdepvar(:,2)-initialdepvar(:,1)-freepar(end)*OCMATHET.initialstatedifference];
else
    resinit=[resinit; ...
        initialdepvar(:,2)-initialdepvar(:,1)];
end
res=[resinit;resasym;resconnec;resequilibrium;resricatti;resper];
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
function [b,infoS,labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATHET OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
sol=evalatmesh(tmesh,y,z);

b=0;
infoS=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
depvara=y(:,leftarcindex);
depvarb=y(:,rightarcindex);
depvaralc=depvara(:,OCMATHET.limitcycleindex==1);
depvarblc=depvarb(:,OCMATHET.limitcycleindex==1);

counter=0;
ctr=1;

for arc=1:OCMATCONT.HE.numarc
    solutionindex=OCMATHET.solutionindex(arc);
    limitcycleindex=OCMATHET.limitcycleindex(arc);
    if solutionindex>1
        relarc=arc-OCMATHET.cumsumnumarc(solutionindex-1);
        transformedtimeshift=OCMATHET.cumsumnumarc(solutionindex-1);%OCMATHET.numarc(solutionindex-1);
    else
        relarc=arc;
        transformedtimeshift=0;
    end
    if limitcycleindex
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(solutionindex))];
    else
        if arc>OCMATHET.cumsumnumarc(solutionindex)
            ctr=ctr+1;
        end
        if OCMATHET.freeendtime(solutionindex)
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
        else
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
        end
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
    reverseflag=0;
    if ~limitcycleindex
        if OCMATHET.solutionindex(arc)==1 && strcmp(OCMATHET.pathtype{1},'u')
            reverseflag=1;
        elseif OCMATHET.solutionindex(arc)==3 && strcmp(OCMATHET.pathtype{2},'u')
            reverseflag=1;
        end
    end
    if ~reverseflag
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
ctr=0;
for ii=1:OCMATHET.hetorder
    %actdepvara=depvara(:,OCMATHET.solutionindex==ii);
    actdepvarb=depvarb(:,OCMATHET.solutionindex==ii);
    limitcycleindex=OCMATHET.limitcycleindex(OCMATHET.solutionindex==ii);
    if ~limitcycleindex(1)
        ctr=ctr+1;
    end
    if ~OCMATHET.freeendtime(ii)
        if ~limitcycleindex(1)
            switch OCMATHET.limitsettype{ctr}
                case 'e'
                    hatx=freepar(OCMATHET.equilibriumcoordinate{ctr});
                case 'l'
                    hatx=depvaralc(:,1);
            end
            violationmat=OCMATCONT.OPTIONS.maxdistance-norm(actdepvarb(:,end)-hatx)<0;
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
    end
end

%----------------------------------------------------------------
function h=plotcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
[s,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
ctr=0;
clf
t=[];
for ii=1:OCMATHET.hetorder
    solutionindex=find(OCMATHET.solutionindex==ii);
    limitcycleindex=OCMATHET.limitcycleindex(solutionindex(1));
    if ~limitcycleindex
        ctr=ctr+1;
        transformedtimeshift=OCMATHET.cumsumnumarc(ctr);
        if OCMATHET.freeendtime(ii)
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{ii}).' freepar(OCMATHET.truncationtimecoordinate(ii))];
        else
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{ii}).' OCMATHET.truncationtime(ii)];
        end
        diffarctime=diff(arctime);
        for arc=1:OCMATHET.numarc(ii)
            t=[t diffarctime(arc)*(s(leftarcindex(solutionindex(arc)):rightarcindex(solutionindex(arc)))-transformedtimeshift)+(arctime(arc)-diffarctime(arc)*(arc-1))];
        end
    end
end
h=OCMATHET.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATHET OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
if isempty(OCMATHET.targetparametervalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(end));
else
    fprintf(1,' Difference to targetvalue: %g\n',OCMATHET.targetparametervalue-modelpar(OCMATHET.targetparameterindex));
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
                if isempty(OCMATHET.targetparametervalue)
                    out=coeff(end);
                else
                    [t,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent);
                    out=OCMATHET.targetparametervalue-modelpar(OCMATHET.targetparameterindex);
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
out.arcinterval=[];
out.arcposition=[OCMATCONT.HE.TIMEDDATA.leftarcindex; ...
    OCMATCONT.HE.TIMEDDATA.rightarcindex];

for solutionindex=1:OCMATHET.hetorder
    limitcycleindex=OCMATHET.limitcycleindex(OCMATHET.solutionindex==solutionindex);
    if limitcycleindex(1)
        arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.periodcoordinate(solutionindex))];
    else
        if OCMATHET.freeendtime(solutionindex)
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' freepar(OCMATHET.truncationtimecoordinate(solutionindex))];
        else
            arctime=[0 freepar(OCMATHET.switchtimecoordinate{solutionindex}).' OCMATHET.truncationtime(solutionindex)];
        end
    end
    out.solverinfo.arcinterval{solutionindex}=arctime;
    out.solverinfo.timehorizon(solutionindex)=arctime(end);
end
out.arcarg=OCMATCONT.HE.arcarg;
out.x0=OCMATHET.initialtime;
%out.modelparameter=OCMATHET.parametervalue;
%out.modelparameter(OCMATHET.parameterindex)=freepar(OCMATHET.parametervaluecoord);
out.modelparameter=modelpar;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.bvpmethod;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.conttype='heteroclinicep2lp';
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATHET.pathtype;
out.solverinfo.multiarccalc=OCBVP.multiarccalc;
out.solverinfo.arcarg=OCMATHET.arcarg;
out.solverinfo.switchtimecoordinate=OCMATHET.switchtimecoordinate;
out.solverinfo.truncationtimecoordinate=OCMATHET.truncationtimecoordinate;
out.solverinfo.equilibriumcoordinate=OCMATHET.equilibriumcoordinate;
out.solverinfo.periodcoordinate=OCMATHET.periodcoordinate;
out.solverinfo.parametervaluecoord=OCMATHET.parametervaluecoord;
out.solverinfo.parameterindex=OCMATHET.parameterindex;
out.solverinfo.solutionindex=OCMATHET.solutionindex;
out.solverinfo.limitcycleindex=OCMATHET.limitcycleindex;
out.solverinfo.monodromy=OCMATCONT.monodromy;

if ~OCMATHET.simple
out.solverinfo.Ycoordinate=OCMATHET.Ycoordinate;
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
global OCMATHET OCBVP

OCBVP.limtcycleindex=OCMATHET.limitcycleindex;

% ------------------------------------------------------

function WorkspaceDone


% ------------------------------------------------------
function [tmesh,y,z,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATHET

%numarc=OCMATCONT.HE.numarc;
%domainddata=OCMATCONT.DOMAINDDATA;
modelpar=OCMATHET.parametervalue;
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
        save([OCMATHET.basicglobalvarfilename '4heteroclinicep2lp'],'MODELINFO')
    end
    save([OCMATHET.basicresultfilename '4heteroclinicep2lp'],'sout','bvpout')
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
