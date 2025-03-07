function out=extremalshoot4ft()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@gradientsolution;
out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
out{11}=@operatorpfrechet;
out{16}=@done;
out{17}=@adapt;
out{18}=@predictextremal;
out{19}=@predictextremaldiff;
out{20}=@workspaceadapt;
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
out{25}=@plotsolution;
out{26}=@saveintermediate;
out{27}=@datapath;
out{30}=@printcontinuation;

function [res,extremal]=operatoreq(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);

ctrl_loc=extremal.v;
depvar0=[extremal.y(:,1);extremal.cst_y(:,1);0];
t=extremal.t;
arcid=extremal.arcid;
switch OCSHOOTCONT.conttype
    case 'time'
        if OCSHOOTCONT.contindex==1
        elseif OCSHOOTCONT.contindex==2
            t=t*freepar(OCSHOOTCONT.continuationindex);
        end
    case 'initialstate'
        depvar0(OCSHOOTCONT.fixinitalstatecoordinate)=OCSHOOTCONT.fixinitialstate+freepar*OCSHOOTCONT.targetvector;
end
depvar0(OCSHOOTCONT.costatecoordinate(OCSHOOTCONT.fixinitialstatecoordinate))=freepar(OCSHOOTCONT.fixinitialstateindex);
[depvarnew,ctrl_locnew,arcidnew]=OCSHOOTCONT.shootingsolver(t,depvar0,ctrl_loc,modelpar,arcid);
res=OCGRADSOL.transversaltycondition(t(end),depvarnew(:,end),modelpar);
res=[res;0];
extremal.y=depvarnew(OCSHOOTCONT.statecoordinate,:);
extremal.cst_y=depvarnew(OCSHOOTCONT.costatecoordinate,:);
extremal.o=depvarnew(OCSHOOTCONT.objectivecoordinate,:);
extremal.v=ctrl_locnew;
extremal.arcid=arcidnew;
[extremal.tangent]=directionalderivative(coeff,extremal,tangent);

function [res,extremal]=operatoreqred(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);

ctrl_loc=extremal.v;
depvar0=[extremal.y(:,1);extremal.cst_y(:,1);0];
t=extremal.t;
arcid=extremal.arcid;
switch OCSHOOTCONT.conttype
    case 'time'
        if OCSHOOTCONT.contindex==1
        elseif OCSHOOTCONT.contindex==2
            t=t*freepar(OCSHOOTCONT.continuationindex);
        end
    case 'initialstate'
        depvar0(OCSHOOTCONT.fixinitalstatecoordinate)=OCSHOOTCONT.fixinitialstate+freepar*OCSHOOTCONT.targetvector;
end
depvar0(OCSHOOTCONT.costatecoordinate(OCSHOOTCONT.fixinitialstatecoordinate))=freepar(OCSHOOTCONT.fixinitialstateindex);
[depvarnew,ctrl_locnew,arcidnew]=OCSHOOTCONT.shootingsolver(t,depvar0,ctrl_loc,modelpar,arcid);
res=OCGRADSOL.transversaltycondition(t(end),depvarnew(:,end),modelpar);
res=[res;0];
extremal.y=depvarnew(OCSHOOTCONT.statecoordinate,:);
extremal.cst_y=depvarnew(OCSHOOTCONT.costatecoordinate,:);
extremal.o=depvarnew(OCSHOOTCONT.objectivecoordinate,:);
extremal.v=ctrl_locnew;
extremal.arcid=arcidnew;

function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;

function J=frechetder(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
J=zeros(OCSHOOTCONT.state_num.concentrated,length(coeff));
depvarT=[extremal.y(:,end);extremal.cst_y(:,end)];
t=extremal.t;

directionalderivativestep=OCSHOOTCONT.OPTIONS.directionalderivativestep;
res0=OCGRADSOL.transversaltycondition(t(end),depvarT,modelpar);
for ii=1:length(coeff)
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;
    res=operatoreqred(coeff1,extremal,tangent);
    J(:,ii)=(res(1:end-1)-res0)/directionalderivativestep;
end

%-------------------------------------------------------------------------
function extremal=predictextremal(coeff,extremal,stepwidth)
global OCSHOOTCONT OCGRADSOL
v0=extremal.v;
Dv=extremal.tangent.v;
offset=0;
for ii=1:length(coeff)
    v0=v0+stepwidth(ii)*Dv(OCSHOOTCONT.loccontrolcoordinate+offset);
    offset=offset+OCSHOOTCONT.loccontrolcoordinate(end);
end
extremal.tangent=[];
extremal.v=v0;
%----------------------------------------------------------------
function h=plotcontinuation(coeff,extremal,tangent,contdata,makemovie)
global OCGRADSOL OCSHOOTCONT
[freepar,modelpar]=drearr(coeff);
t=extremal.t;
switch OCSHOOTCONT.conttype
    case 'time'
        if OCSHOOTCONT.contindex==1
        elseif OCSHOOTCONT.contindex==2
            t=t*freepar(OCSHOOTCONT.continuationindex);
        end
    case 'initialstate'
        depvar0(OCSHOOTCONT.fixinitalstatecoordinate)=OCSHOOTCONT.fixinitialstate+freepar*OCSHOOTCONT.targetvector;
end
extremal.t=t;
h=OCGRADSOL.plotcontinuation(freepar,extremal,tangent,modelpar,contdata,OCSHOOTCONT.conttype,makemovie);

%----------------------------------------------------------------
function h=plotsolution(coeff,extremal,tangent)
global OCGRADSOL OCSHOOTCONT
[freepar,modelpar]=drearr(coeff);
h=OCGRADSOL.plotsolution(coeff,extremal,tangent,modelpar);


%----------------------------------------------------------------
function idx=printcontinuation(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL
idx=[];
if isempty(coeff)
    return
end
if isempty(OCSHOOTCONT.targetvalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCSHOOTCONT.continuationindex));
else
    out=targetvaluefunc(1,coeff,extremal,tangent);
    fprintf(1,' Distance from Target value: %g\n',out);
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCSHOOTCONT OCGRADSOL

varargout{2}=nan;
% all done succesfully
varargout{1}=0;
%-------------------------------------------------------------
function [out, failed]=testfunc(id,coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL

out(1)=0;
failed=[];
%----------------------------------------------------------------
function [out, failed]=targetvaluefunc(id,coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL

failed=[];
for ii=id
    switch ii
        case 1
            if isempty(OCSHOOTCONT.targetvalue)
                out=[];
            else
                switch OCSHOOTCONT.conttype
                    case 'time'
                        freepar=drearr(coeff);
                        out=OCSHOOTCONT.targetvalue-freepar(end);
                    case 'initialstate'
                        freepar=drearr(coeff);

                        out=1-freepar;
                    case 'parameter'
                        freepar=drearr(coeff);

                        out=OCSHOOTCONT.targetvalue-freepar;
                end
            end

        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function b=testadmissibility(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL 
[freepar,modelpar]=drearr(coeff);
switch OCSHOOTCONT.conttype
    case 'time'
        if OCSHOOTCONT.contindex==1
        elseif OCSHOOTCONT.contindex==2
            extremal.t=extremal.t*freepar;
            extremal.timehorizon=freepar;
        end
end
val=OCSHOOTCONT.admissibility(extremal.t,extremal.y,extremal.v,modelpar);
b=any(val(:)<0);

%----------------------------------------------------------------
function out=formatsolution(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL 
[freepar,modelpar]=drearr(coeff);
t=extremal.t;
switch OCSHOOTCONT.conttype
    case 'time'
        if OCSHOOTCONT.contindex==1
        elseif OCSHOOTCONT.contindex==2
            t=t*freepar(OCSHOOTCONT.continuationindex);
        end
    case 'initialstate'
        depvar0(OCSHOOTCONT.fixinitalstatecoordinate)=OCSHOOTCONT.fixinitialstate+freepar*OCSHOOTCONT.targetvector;
end
extremal.t=t;
out.extremal=extremal;
out.parameters=freepar;
out.extremal.modelparameter=modelpar;
out.extremal.modelname=OCSHOOTCONT.modelname;
out.extremal.type=OCSHOOTCONT.octype;
out.extremal.contclass=OCSHOOTCONT.problem_func;
out.extremal.conttype=OCSHOOTCONT.conttype;
switch OCSHOOTCONT.conttype
    case 'time'
    case 'initialstate'
    case 'parameter'
        out.extremal.freeparameterindex=OCSHOOTCONT.contindex;
end
%---------------------------------------------------------
function done
WorkspaceDone;
% ---------------------------------------------------------

function WorkspaceInit(coeff,extremal,tangent)
global OCSHOOTCONT OCGRADSOL

% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------

function [freepar,modelpar]=drearr(coeff)
global OCSHOOTCONT OCGRADSOL

modelpar=OCSHOOTCONT.modelparameter;
freepar=coeff;
switch  OCSHOOTCONT.conttype
    case 'parameter'
        modelpar(OCSHOOTCONT.freeparameterindex)=freepar(OCSHOOTCONT.continuationindex);%OCSHOOTCONT.startingparametervalue+freepar*OCSHOOTCONT.targetvector;
end

%-----------------------------------------------------------------
function failed=saveintermediate(sout,gradout,contnum)
global OCSHOOTCONT OCGRADSOL
failed=0;
MODELINFO.OCSHOOTCONT=OCSHOOTCONT;
MODELINFO.OCGRADSOL=OCGRADSOL;
try
    if contnum==1
        save([OCGRADSOL.basicglobalvarfilename '4extremalshoot4ft'],'MODELINFO')
    end
    save([OCGRADSOL.basicresultfilename '4extremalshoot4ft'],'sout','gradout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCGRADSOL

pathname=OCGRADSOL.datapath();

% helper function for directional derivative
function dX=directionalderivative(coeff,extremal,tangent)
global OCSHOOTCONT

directionalderivativestep=OCSHOOTCONT.OPTIONS.directionalderivativestep;
y=extremal.y;
v=extremal.v;
cst_y=extremal.cst_y;
dX.y=[];
dX.v=[];
dX.cst_y=[];
for ii=1:length(coeff)
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;
    [dum,extremal1]=operatoreqred(coeff1,extremal,tangent);
    dx.y=(extremal1.y-y)/directionalderivativestep;
    dx.v=(extremal1.v-v)/directionalderivativestep;
    dx.cst_y=(extremal1.cst_y-cst_y)/directionalderivativestep;
    dX.y=[dX.y;dx.y];
    dX.cst_y=[dX.cst_y;dx.cst_y];
    dX.v=[dX.v;dx.v];
end

