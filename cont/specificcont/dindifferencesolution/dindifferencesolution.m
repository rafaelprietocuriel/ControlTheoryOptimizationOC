function out=dindifferencesolution()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@odiff;
out{4}{2}=@bc;
out{4}{3}=@ic;
out{5}{1}=@odejac;
out{5}{2}=@bcjac;
out{5}{3}=@icjac;
out{6}=@findarcposition;
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
out{30}=@printcontinuation;

function res=operatoreq(tmesh,coeff,tangent,odefun,bcfun,icfun)
global OCMATCONT
[t,y,freepar,modelpar]=drearr(tmesh,coeff);
res=calc_DMRHS(t,y,freepar,modelpar,odefun,bcfun,icfun);
res(OCMATCONT.HE.numdvariables,1)=0;

function J=frechetder(tmesh,coeff,tangent,odefun,bcfun,icfun,odejacfun,bcjacfun,icjacfun)

[t,y,freepar,modelpar]=drearr(tmesh,coeff);
J=calc_DMRHSJac(t,y,freepar,modelpar,odefun,bcfun,icfun,odejacfun,bcjacfun,icjacfun);

% for testing
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun,icfun},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% J=full(J);
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end
% J=Jnum;

function opt=options(opt)

function varargout=probleminit(varargin)
tmesh=varargin{1};
coeff=varargin{2};
tangent=varargin{3};
WorkspaceInit(tmesh,coeff,tangent);

% all done succesfully
varargout{1} = 0;


%-------------------------------------------------------------------------
function [out,failed]=findarcposition(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF OCBVP
[t,y,freepar,modelpar]=drearr(tmesh,coeff);
out=0;
failed=0;
newarcarg=[];
newarcposition=[];
for ii=1:OCMATINDIF.indifferenceorder
    idx=OCMATINDIF.pathcoord(1,ii):OCMATINDIF.pathcoord(2,ii);
    [arcposition arcarg]=OCMATINDIF.findarcposition(t(idx),y(:,idx),modelpar);
    if isempty(arcposition)
        return
    end
    newarcarg=[newarcarg arcarg];
    newnumarc(ii)=length(arcarg);
    if ii==1
        arcposabs=0;
    else
        arcposabs=arcposition(2,end);
    end
    %arcposition(1,2:end)=arcposition(1,2:end)-1;
    arcposition=arcposabs+arcposition;
    %arcposition(2,end)=arcposition(2,end)+1;
    newarcposition=[newarcposition arcposition];
end
OCMATCONT.arcarg=newarcarg;
OCBVP.Lidx=newarcposition(1,:);
OCBVP.Ridx=newarcposition(2,:);
OCBVP.numarc=newnumarc;
OCBVP.Nint=OCBVP.Ridx-OCBVP.Lidx;
OCBVP.arcposition=newarcposition;
OCMATCONT.HE.arcindex=arcarg2arcindex(newarcarg);
OCMATCONT.HE.arcarg=newarcarg;

%-------------------------------------------------------------------------
function Dx=odiff(t,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
arcarg=OCMATCONT.HE.arcarg(arc);
Dx=OCMATINDIF.canonicalsystemmap(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
% Remark: for a nonautonomous problem the Jacobian at the switching times
% has to be adapted; it has to include derivatives of the form dtds*f_t
function [J Jpar]=odejac(t,depvar,arc,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
arcarg=OCMATCONT.HE.arcarg(arc);
J=OCMATINDIF.canonicalsystemmapjacobian(t,depvar,modelpar,arcarg);
Jpar=zeros(OCBVP.neqn,OCBVP.npar);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP
resinit=[];
resasym=[];

%restarget=OCMATINDIF.bcinitial(depvara(:,1),OCMATINDIF.targetcoordinate,OCMATINDIF.startvalue+freepar(end)*OCMATINDIF.continuationvector);
for order=1:OCMATINDIF.indifferenceorder
    if order==1
        if ~isempty(OCMATINDIF.freevector)
            initialstate=OCMATINDIF.startvalue;
            for jj=1:length(OCMATINDIF.freevectorcoord)
                initialstate=initialstate+freepar(OCMATINDIF.freevectorcoord(jj))*OCMATINDIF.freevector(:,jj);
            end
            restarget=depvara(OCMATINDIF.statecoordinate,1)-initialstate;
        end
    end
    if order<OCMATINDIF.indifferenceorder
        resinit=[resinit; ...
            depvara(OCMATINDIF.statecoordinate,OCMATINDIF.arcargcoord(1,order+1))-depvara(OCMATINDIF.statecoordinate,OCMATINDIF.arcargcoord(1,order))];
    end
    resasym=[resasym;OCMATINDIF.bcasymptotic(depvarb(:,OCMATINDIF.arcargcoord(2,order)),OCMATINDIF.asymptoticmatrix{order},OCMATINDIF.saddlepoint{order})];
end
res=[restarget;resinit;resasym];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATINDIF OCMATCONT OCBVP


Jpar=zeros(OCBVP.nBCs,OCMATCONT.HE.numparameter);
Ja=zeros(OCBVP.nBCs,sum(OCBVP.numarc)*OCBVP.nummap);
Jb=Ja;
colcounter=0;
rowcounter=0;
colidx_start=colcounter+1;
colcounter=colcounter+OCBVP.nummap;
[Japart Jbpart Jpartmp]=OCMATINDIF.bcjacobianinitial(depvara,freepar,OCMATCONT.HE.arcarg(1),OCMATINDIF.targetcoordinate,OCMATINDIF.continuationvector);
Jpar(1:OCMATCONT.HE.numinitialcondition,:)=Jpartmp;
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numinitialcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;
colcounterb=0;
colcounter=0;
for order=1:OCMATINDIF.indifferenceorder
    if order<OCMATINDIF.indifferenceorder
        rowcounterstart=rowcounter+1;
        rowcounter=rowcounter+OCMATINDIF.statenum;
        eyen=eye(OCMATINDIF.statenum);
        colidx_start=colcounter+1;
        colcounter=colcounter+OCMATINDIF.statenum;
        Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=-eyen;
        colcounter=colcounter+OCMATINDIF.statenum;
        colidx_start=colcounter+(OCBVP.numarc(order)-1)*OCBVP.nummap+1;
        colcounter=colcounter+(OCBVP.numarc(order)-1)*OCBVP.nummap+OCMATINDIF.statenum;
        Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=eyen;
        colcounter=colcounter+OCMATINDIF.statenum;
    end
    [Japart Jbpart]=OCMATINDIF.bcjacobianasymptotic(depvarb,OCMATCONT.HE.arcarg(OCBVP.numarc(order)),OCMATINDIF.asymptoticmatrix{order}',OCMATINDIF.saddlepoint{order});
    if OCBVP.numarc(order)>1
        colcounterb=colcounterb+(OCBVP.numarc(order)-1)*OCBVP.nummap;
    end
    rowcounterstart=rowcounter+1;
    rowcounter=rowcounter+OCMATINDIF.saddlepointcodimension(order);
    colidx_start=colcounterb+1;
    colcounterb=colcounterb+OCBVP.nummap;
    Jb(rowcounterstart:rowcounter,colidx_start:colcounterb)=Jbpart;
end

%-------------------------------------------------------------------------
function res=ic(t,depvar,freepar,modelpar)
global OCMATINDIF

for ii=1:OCMATINDIF.indifferenceorder
    arcpos=OCMATINDIF.arcposition(:,OCMATINDIF.arcargcoord(1,ii):OCMATINDIF.arcargcoord(2,ii));
    arcarg=OCMATINDIF.arcarg(OCMATINDIF.arcargcoord(1,ii):OCMATINDIF.arcargcoord(2,ii));
    totalo=[];
    for jj=1:OCMATINDIF.arcnum(ii)
        arcp=arcpos(1,jj):arcpos(2,jj);
        o=OCMATINDIF.objectivefunction(t(arcp(2:end)),depvar(:,arcp),modelpar,arcarg(jj));
        totalo(arcp(1:end-1))=o;
    end
    O(ii)=sum(totalo);
end
res=diff(O);

%-------------------------------------------------------------------------
function [J Jpar]=icjac(t,depvar,freepar,modelpar)
global OCMATINDIF
Jpar=zeros(1,length(freepar));

sign=1;
for ii=1:OCMATINDIF.indifferenceorder
    sign=-sign;
    arcpos=OCMATINDIF.arcposition(:,OCMATINDIF.arcargcoord(1,ii):OCMATINDIF.arcargcoord(2,ii));
    arcarg=OCMATINDIF.arcarg(OCMATINDIF.arcargcoord(1,ii):OCMATINDIF.arcargcoord(2,ii));
    for jj=1:OCMATINDIF.arcnum(ii)
        arcp=arcpos(1,jj):arcpos(2,jj);
        [Jxtmp Jltmp]=OCMATINDIF.objectivefunctionjacobian(t(arcp(2:end)),depvar(:,arcp),modelpar,arcarg(jj));
        if sign==-1
            Jx(:,arcp(1:end-1))=-Jxtmp;
            Jl(:,arcp(1:end-1))=-Jltmp;
        else
            Jx(:,arcp(1:end-1))=Jxtmp;
            Jl(:,arcp(1:end-1))=Jltmp;
        end
    end
end
zerovec=zeros(size(Jx,1),1);
J=[Jx zerovec;zerovec Jl];
J=J(:).';
%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATINDIF OCBVP

[t,y,freepar,modelpar]=drearr(tmesh,coeff);

b=0;
infoS=[];
leftarcindex=OCBVP.Lidx;
rightarcindex=OCBVP.Ridx;
counter=0;
for arc=1:OCBVP.numarc
    arcindex=OCMATCONT.HE.arcindex(arc);
    arcarg=OCMATCONT.HE.arcarg(arc);
    t=tmesh(leftarcindex(arc):rightarcindex(arc));
    [constr labelS]=OCMATINDIF.testadmissibility(t,y(OCBVP.cols,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
global OCMATINDIF OCMATCONT
[t,y,freepar,modelpar]=drearr(tmesh,coeff,tangent);

h=OCMATINDIF.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
drawnow
figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
%global OCMATINDIF OCMATCONT
idx=[];
if isempty(coeff)
    return
end
[t,y,freepar,modelpar]=drearr(tmesh,coeff,tangent);
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
global OCMATCONT

failed=[];
for ii=id
    switch ii
        case 1
            out=OCMATCONT.OPTIONS.totalrelativedistance-coeff(OCMATCONT.HE.contparametercoord);
        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF OCBVP
[t,y,freepar,modelpar]=drearr(tmesh,coeff,tangent);

out.x0=tmesh(1);
out.y0=y(:,1);
out.x=tmesh(2:end);
out.y=y(:,2:end);

out.arcarg=OCMATCONT.HE.arcarg;
out.arcposition=OCBVP.arcposition;
out.arcposition(2,:)=out.arcposition(2,:)-1;
out.timehorizon=repmat(inf,1,OCMATINDIF.indifferenceorder);
out.modelparameter=OCMATINDIF.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.newtonsolver;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.arcnum=OCMATINDIF.arcnum;
out.solverinfo.pathcoord=OCMATINDIF.pathcoord;
out.solverinfo.arcargcoord=OCMATINDIF.arcargcoord;
out.solverinfo.parameters=freepar;
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATINDIF.pathtype;

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATINDIF
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent);
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
function [tmesh,y,freepar,modelpar]=drearr(tmesh,coeff,tangent)
global OCMATCONT OCMATINDIF

modelpar=OCMATINDIF.parametervalue;
y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
z=[];
freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values

%-----------------------------------------------------------------
function failed=saveintermediate(sout,bvpout,contnum)
global OCMATCONT OCMATINDIF OCBVP
failed=0;
MODELINFO.OCMATCONT=OCMATCONT;
MODELINFO.OCMATINDIF=OCMATINDIF;
MODELINFO.OCBVP=OCBVP;
try
    if contnum==1
        save([OCMATINDIF.basicglobalvarfilename '4dindifferencesolution'],'MODELINFO')
    end
    save([OCMATINDIF.basicresultfilename '4dindifferencesolution'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATINDIF

pathname=OCMATINDIF.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;
