function out=dextremal2fp()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@options;
out{4}{1}=@odiff;
out{4}{2}=@bc;
out{4}{3}=[];
out{5}{1}=@odiffjac;
out{5}{2}=@bcjac;
out{5}{3}=[];
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
% Jnum=numjaccsd(@operatoreq,{tmesh,coeff,tangent,odefun,bcfun},numel(coeff),numJacOpt);
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
global OCMATCONT OCMATAE OCBVP
out=0;
failed=0;
[t,y,freepar,modelpar]=drearr(tmesh,coeff);
[arcposition arcarg]=OCMATAE.findarcposition(t,y,modelpar);
if isempty(arcposition)
    return
end
OCMATCONT.arcarg=arcarg;
OCBVP.Lidx=arcposition(1,:);
OCBVP.Ridx=arcposition(2,:);
OCBVP.numarc=numel(arcarg);
if OCBVP.numarc>1 
    OCBVP.multipointbvp=true;
else
    OCBVP.multipointbvp=false;
end
OCBVP.Nint=OCBVP.Ridx-OCBVP.Lidx;
OCBVP.arcposition=arcposition;
OCMATCONT.HE.arcindex=arcarg2arcindex(arcarg);
OCMATCONT.HE.arcarg=arcarg;

%-------------------------------------------------------------------------
function [F J]=operatorpfrechet(tmesh,coeff,tangent,varargin)
F=operatoreq(tmesh,coeff,tangent,varargin{:});
J=frechetder(tmesh,coeff,tangent,varargin{:});

%-------------------------------------------------------------------------
function res=odiff(t,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT
arcarg=OCMATCONT.HE.arcarg(arc);
res=OCMATAE.canonicalsystemmap(t,depvar,modelpar,arcarg);

%-------------------------------------------------------------------------
function [J Jpar]=odiffjac(t,depvar,arc,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
arcarg=OCMATCONT.HE.arcarg(arc);
J=OCMATAE.canonicalsystemmapjacobian(t,depvar,modelpar,arcarg);
Jpar=zeros(OCBVP.neqn,OCBVP.npar);

%-------------------------------------------------------------------------
function res=bc(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP
if strcmp(OCMATAE.pathtype,'u')
    tmp=depvarb;
    depvarb=depvara(:,end:-1:1);
    depvara=tmp(:,end:-1:1);
end

initialstate=OCMATAE.startvalue+freepar(end)*OCMATAE.continuationvector;
resinit=OCMATAE.bcinitial(depvara,OCMATAE.targetcoordinate,initialstate,modelpar,OCMATCONT.HE.arcarg(1));
resasym=OCMATAE.bcasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint);
res=[resinit;resasym];

%-------------------------------------------------------------------------
function [Ja Jb Jpar]=bcjac(depvara,depvarb,freepar,modelpar)
global OCMATAE OCMATCONT OCBVP

Jpar=zeros(OCBVP.nBCs,OCMATCONT.HE.numparameter);
Ja=zeros(OCBVP.nBCs);
Jb=Ja;
if strcmp(OCMATAE.pathtype,'u')
    tmp=depvarb;
    depvarb=depvara(:,end:-1:1);
    depvara=tmp(:,end:-1:1);
end
colcounter=0;
rowcounter=0;
colidx_start=colcounter+1;
colcounter=colcounter+OCBVP.nummap;
[Japart Jbpart Jpartmp]=OCMATAE.bcjacobianinitial(depvara,freepar,OCMATCONT.HE.arcarg(1),OCMATAE.targetcoordinate,OCMATAE.continuationvector);
Jpar(1:OCMATCONT.HE.numinitialcondition,:)=Jpartmp;
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numinitialcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;
if OCBVP.numarc>1
     colidx_start=colcounter+(OCBVP.numarc-2)*OCBVP.nummap+1;
     colcounter=colcounter+(OCBVP.numarc-1)*OCBVP.nummap;
end
[Japart Jbpart]=OCMATAE.bcjacobianasymptotic(depvarb,OCMATCONT.HE.arcarg(OCBVP.numarc),OCMATAE.asymptoticmatrix',OCMATAE.saddlepoint);
rowcounterstart=rowcounter+1;
rowcounter=rowcounter+OCMATCONT.HE.numendcondition;
Ja(rowcounterstart:rowcounter,colidx_start:colcounter)=Japart;
Jb(rowcounterstart:rowcounter,colidx_start:colcounter)=Jbpart;
if strcmp(OCMATAE.pathtype,'u')
    I=reshape(1:colcounter,OCBVP.nummap,[])';
    I=I(end:-1:1,:)';
    I=reshape(I,1,[]);
    tmp=Jb;
    Jb=Ja(:,I);
    Ja=tmp(:,I);
end

%----------------------------------------------------------------
function [b infoS labelS]=testadmissibility(tmesh,coeff,tangent,varargin)
global OCMATCONT OCMATAE OCBVP

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
    [constr labelS]=OCMATAE.testadmissibility(t,y(OCBVP.cols,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg);
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
[t,y,freepar,modelpar]=drearr(tmesh,coeff);

h=OCMATAE.plotcontinuation(t,y,modelpar,OCMATCONT.HE.arcarg,freepar,tangent);
%set(gcf,'WindowStyle','docked')
% drawnow
% figure(gcf)

%----------------------------------------------------------------
function idx=printcontinuation(tmesh,coeff,tangent)
global OCMATCONT
idx=[];
if isempty(coeff)
    return
end
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
global OCMATCONT OCMATAE OCBVP
[t,y,freepar]=drearr(tmesh,coeff);

switch OCMATAE.pathtype
    case 's'
        out.x0=tmesh(1);
        out.y0=y(:,1);
        out.x=tmesh(2:end);
        out.y=y(:,2:end);
    case 'u'
        out.x0=tmesh(end);
        out.y0=y(:,end);
        out.x=tmesh(1:end-1);
        out.y=y(:,1:end-1);
end

out.arcarg=OCMATCONT.HE.arcarg;
out.arcposition=OCBVP.arcposition;
out.arcposition(2,:)=out.arcposition(2,:)-1;
out.timehorizon=inf;
out.modelparameter=OCMATAE.parametervalue;
out.modelname=OCMATCONT.modelname;

out.solver=OCMATCONT.newtonsolver;
out.solverinfo.coeff=coeff;
out.solverinfo.tangent=tangent;
out.solverinfo.tmesh=tmesh;
out.solverinfo.parameters=freepar;
out.solverinfo.continuationparameter=coeff(OCMATCONT.HE.contparametercoord);
out.solverinfo.pathtype=OCMATAE.pathtype;

%---------------------------------------------------------------------
function [failed,s]=process(id,tmesh,coeff,tangent,s)
global OCMATCONT OCMATAE
fprintf('\n label=%s\n Continuation parameter=%g\n', s.label,coeff(OCMATCONT.HE.contparametercoord));
switch id
    case 1 % LP
        J=frechetder(tmesh,coeff,tangent,OCMATCONT.odiff,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac);
        [Q,R,E]=qr(full(J)); % see Govaerts et al 2005 (p. 242)
        R(end,end)=0;R(end,end-1)=0;
        p=E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
        OCMATAE.PD_phi=p'/norm(p);
        OCMATAE.PD_psi=Q(:,end);
        s.data.phi=OCMATAE.PD_phi(:);
        s.data.laecoefficient=[];%nf_LAE(tmesh,coeff); % quadratic coefficient of center manifold
        s.data.sol=formatsolution(tmesh,coeff,tangent);
        s.msg =sprintf('Limit asymptotic dextremal');
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
function [tmesh,y,freepar,modelpar]=drearr(tmesh,coeff)
global OCMATCONT OCMATAE

modelpar=OCMATAE.parametervalue;
y=coeff(OCMATCONT.HE.DDATA.meshvalcoord);
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
        save([OCMATAE.basicglobalvarfilename '4dextremal2fp'],'MODELINFO')
    end
    save([OCMATAE.basicresultfilename '4dextremal2fp'],'sout','bvpout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCMATAE

pathname=OCMATAE.datapath();

%-----------------------------------------------------------------
function [flag,tmesh,coeff,tangent] = adapt(tmesh,coeff,tangent)
flag=0;
